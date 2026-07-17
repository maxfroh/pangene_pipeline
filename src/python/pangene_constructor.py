#!/usr/bin/env python3
import bisect
import csv
import gzip
import logging
import queue
import re
import shutil
import sqlite3 as sql
import string
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import (ProcessPoolExecutor, ThreadPoolExecutor,
                                as_completed)
from concurrent.futures.process import BrokenProcessPool
from itertools import combinations, product
from pathlib import Path
from typing import Any, Generator

import gffutils
import igraph as ig
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pyfaidx import Fasta
from tqdm import tqdm

from .param_manager import ParamManager
from .utils import (build_logger, concat_files, copy_file_quiet, execute,
                    execute_quiet, get_name_ext_and_is_gzip, gunzip_file_quiet,
                    strip_filename)


class PangeneConstructor:
    """A factory class which creates a pangene reference."""

    def __init__(
        self,
        reference: str,
        pangenes_dir: Path,
        pangene_info: dict[str, str | float],
        pm: ParamManager,
    ):
        """
        Constructs a pangene reference.

        :param reference: The name of the pangene reference.
        :type reference: str
        :param pangenes_dir: The directory to store the constructed pangene.
        :type pangenes_dir: Path
        :param pangene_info: The dictionary containing all pangene-specific information provided by the user.
        :type pangene_info: dict[str, str | float]
        :param pm: The `ParamManager` object with all pipeline parameters.
        :type pm: ParamManager
        """
        self.pm = pm
        self.reference = reference
        self.constructed = False
        self.logger = build_logger(f"{str(self)}")

        self.pangene_dir = pangenes_dir.resolve() / reference
        self.pangene_dir.mkdir(exist_ok=True, parents=True)
        self.tmp_dir = self.pangene_dir / "tmp"
        self.tmp_dir.mkdir(exist_ok=True, parents=True)

        try:
            self.grp_file = pangene_info.get("grp_file", None)
            if self.grp_file:
                self.grp_file = Path(self.grp_file)
            self.pangene_fastas_dir = Path(pangene_info["pangene_fastas_dir"])
            self.add_k_vals_pl = Path(pangene_info["add_k_vals_pl"])
            self.redunancy_thresh = pangene_info["redundancy_thresh"]
            self.max_diff_r = pangene_info.get("max_diff_r", 0.001)
            ks_threshold_pairs = pangene_info.get("ks_threshold_pairs", [])
            self.ks_threshold_pairs = [
                "_".join(sorted(pair)) for pair in ks_threshold_pairs
            ]
        except KeyError as ke:
            key = str(ke).strip("'").strip('"')
            self.logger.error(
                f"Required [input] key {key} is missing from the configuration!"
            )
            raise ke
        self.annotation_file = self.pangene_dir / "annotation.map"
        self.cds_fasta = self.pangene_dir / f"{reference}_cds.fa.gz"

        self.constituents = [c.name for c in self.pangene_fastas_dir.iterdir()]

        if self.redunancy_thresh < 0.75 or self.redunancy_thresh > 1.0:
            thresh = max(min(0.75, self.redunancy_thresh), 1.0)
            self.logger.info(
                f"Pangene redundancy threshold for {self.reference} is out of bounds! New threshold is {thresh} (was {self.redunancy_thresh})"
            )
            thresh = self.redunancy_thresh

        # Get word_size; see https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#user-content-CDHITEST for source of numbers
        threshholds = [0.8, 0.85, 0.88, 0.9, 0.925, 0.95, 0.975, 1.01]
        word_sizes = [4, 5, 6, 7, 8, 9, 10, 11]
        self.word_size = word_sizes[bisect.bisect(threshholds, self.redunancy_thresh)]

        self.target_pep_fastas_dir = self.tmp_dir / "pep_fastas"
        self.target_pep_fastas_dir.mkdir(exist_ok=True, parents=True)
        self.target_cds_fastas_dir = self.tmp_dir / "cds_fastas"
        self.target_cds_fastas_dir.mkdir(exist_ok=True, parents=True)
        self.original_gffs = list(self.pangene_fastas_dir.glob("*/*.gff3*"))
        self.gffs_dir = self.tmp_dir / "gffs"
        self.gffs_dir.mkdir(exist_ok=True, parents=True)
        self.mcscanx_dir = self.tmp_dir / "mcscanx_dir"
        self.mcscanx_dir.mkdir(exist_ok=True, parents=True)
        self.blastps_dir = self.tmp_dir / "blastps"
        self.blastps_dir.mkdir(exist_ok=True, parents=True)
        self.pangene_db = self.tmp_dir / f"{self.reference}.db"

        self.plots_dir = self.pangene_dir / "plots"
        self.plots_dir.mkdir(exist_ok=True, parents=True)

        codes = [a + b for a, b in product(string.ascii_lowercase, repeat=2)]
        species = [d.name for d in self.pangene_fastas_dir.iterdir() if d.is_dir()]
        self.species_to_code_map = {species[i]: codes[i] for i in range(len(species))}
        self.code_pairs = list(combinations(self.species_to_code_map.values(), 2))

    def get_reference_info(self) -> tuple[Path, Path]:
        if not self.constructed:
            self.logger.info(f"Pangene {str(self)} not constructed! Constructing now.")
            self.construct_pangene()
        return (self.annotation_file, self.cds_fasta)

    def construct_pangene(self):
        """
        if self.grp_file is None:
            # Get necessary FASTAs (peptide/amino acid and CDS)
            self.gather_fastas()
            # Run OrthoFinder
            og_of_file = self.get_orthologous_groups_with_orthofinder()
            # Generate MCScanX-compatible GFF files
            self.prepare_gffs_for_mcscanx()
            # Generate MCScanX blastp information
            self.blastp_alignment_and_filtering()
            # Run MCScanX and add Ks/Ka microsynteny information
            self.run_mcscanx_and_get_ks()
            # Filter MCScanX synteny blocks by a calculated Ks value
            self.calculate_ks_threshold_and_filter()
            # Use MCScanX microsynteny information to subdivide OrthoFinder groups
            new_og_file = self.create_new_orthogroups_with_synteny_info(og_of_file)
            # Turn these subdivisions into new orthologous groups
            sorted_og_df_file = self.rename_orthogroups_and_remove_empty(new_og_file)
            # Create a melted orthologous group file
            self.grp_file = self.melt_orthogroups(sorted_og_df_file)
        # Prune redunant genes from the orthogroups
        self.prune_redundancies()
        # Plot the pruning results
        self.plot_count_difference()
        """
        # Mark that the pangene has been built
        self.constructed = True
        # Housekeeping to save space
        # if self.tmp_dir.exists():
        #     shutil.rmtree(self.tmp_dir)

    def gather_fastas(self):
        """
        Collect both peptide (amino acid) and CDS (coding sequence) FASTA files for use in pangene construction.
        """
        self.logger.info("Collecting peptide and CDS FASTA files.")

        fasta_dirs = [
            d
            for d, _, _ in self.pangene_fastas_dir.walk()
            if d != self.pangene_fastas_dir
        ]
        to_gunzip = []
        to_copy = []
        # For each genome subdirectory
        for fasta_dir in fasta_dirs:
            # For each peptide and CDS FASTA in the subdirectory (should only be one of each)
            for fasta_type, target_fasta_dir in [
                ("pep", self.target_pep_fastas_dir),
                ("cds", self.target_cds_fastas_dir),
            ]:
                genome = fasta_dir.name
                new_file = target_fasta_dir / f"{genome}.fa"
                fasta_file = list(fasta_dir.glob(f"{genome}_{fasta_type}.fa"))
                gunzipped = False
                if fasta_file:
                    gunzipped = True
                else:
                    fasta_file = list(fasta_dir.glob(f"{genome}_{fasta_type}.fa.gz"))
                    if not fasta_file:
                        continue
                fasta_file = fasta_file[0]
                if not gunzipped:
                    to_gunzip.append({"source": fasta_file, "dest": new_file})
                else:
                    to_copy.append({"source": fasta_file, "dest": new_file})

        # Copy files to new directories in parallel
        # TODO: is self.pm.p the best?
        with ProcessPoolExecutor(max_workers=self.pm.p) as executor:
            futures = [
                executor.submit(gunzip_file_quiet, f["source"], f["dest"])
                for f in to_gunzip
            ]
            for future in as_completed(futures):
                result = future.result()
                if result != None:
                    _, source, dest = result
                    self.logger.error(f"Unable to unzip {source} to {dest}.")

        with ProcessPoolExecutor(max_workers=self.pm.p) as executor:
            futures = [
                executor.submit(copy_file_quiet, f["source"], f["dest"])
                for f in to_copy
            ]
            for future in as_completed(futures):
                result = future.result()
                if result != None:
                    _, source, dest = result
                    self.logger.error(f"Unable to copy {source} to {dest}.")

        self.logger.info("Peptide and CDS FASTA files collected successfully!")

    def get_orthologous_groups_with_orthofinder(self) -> Path:
        """
        Run OrthoFinder and write all found orthologous groups to a TSV file.
        :return: The TSV file containing OrthoFinder's results.
        :rtype: Path
        """
        self.logger.info(
            f"Finding orthologous groups for {self.reference} using OrthoFinder."
        )
        fastas_dir = (
            self.target_pep_fastas_dir
        )  # Directory containing FASTA format proteomes to use

        ortho_results_dir = self.tmp_dir / "ortho_results"  # Output directory
        if ortho_results_dir.exists():
            shutil.rmtree(ortho_results_dir)

        # Do not let the number of threads exceed twice the number of FASTA files to prevent race conditions, thrashing, etc.
        max_threads = len(list(fastas_dir.iterdir())) * 2
        # Leave a few cores/threads available if possible to not overwhelm the system
        t = min(max(1, self.pm.p - 4), max_threads)
        # Follow guidance of OrthoFinder manual, a = t/8
        a = max(1, t // 8)

        cmds = [
            "orthofinder",
            "-y",
            "-f",
            fastas_dir,
            "-o",
            ortho_results_dir,
            "-n",
            self.reference,
            "-og",
            "-M",
            "dendroblast",
            "-t",
            t,
            "-a",
            a,
        ]
        execute(
            cmds,
            f"Running OrthoFinder. May take a long time.",
        )

        # Get both assigned and unassigned orthologous groups and combine them
        orthogroups_dir = (
            ortho_results_dir / f"Results_{self.reference}" / "Orthogroups"
        )
        orthogroups = pd.read_csv(orthogroups_dir / "Orthogroups.tsv", sep="\t")
        orthogroups_unassigned = pd.read_csv(
            orthogroups_dir / "Orthogroups_UnassignedGenes.tsv", sep="\t"
        )
        orthogroups_combined = pd.concat([orthogroups, orthogroups_unassigned], axis=0)
        orthogroups_combined = orthogroups_combined.rename(
            columns={"Orthogroup": "OGID"}
        )
        # Genes will be semicolon-delimited, not comma-delimited
        orthogroups_combined = orthogroups_combined[
            [col for col in orthogroups_combined.columns if col != "OGID"]
        ].replace(", ", ";", regex=True)
        of_og_file = self.tmp_dir / "of_orthogroups.tsv"
        orthogroups_combined.to_csv(of_og_file, sep="\t", index=False)

        self.logger.info("OrthoFinder finished successfully!")
        return of_og_file

    @staticmethod
    def _prepare_individual_gff_for_mcscanx(
        gff_file: Path, species_prefix: str, out_dir: Path, gffs_dir: Path
    ):
        """
        Create an MCScanX-compatible GFF file. This is a tab-separated file with four columns.
        * The first column is a two-letter and one-to-two digit string indicating the location of the gene, where
        the two letters correspond to the species (for the pipeline, an assigned code mapped to each species, e.g. `aa1`) and the digits
        correspond to the chromosome number.
        * The second column is the gene identifier.
        * The third column is the start location of the gene.
        * The fourth column is the end location of the gene.

        See https://github.com/wyp1125/MCScanX#mcscanx-1 for more info.

        :param gff_file: The GFF or GFF3 file to use.
        :type gff_file: Path
        :param species_prefix: The code of the species (e.g., `aa`).
        :type species_prefix: str
        :param out_dir: The directory to store the temporary GFF databases.
        :type out_dir: Path
        :param gffs_dir: The directory in which to output the MCScanX-compatible GFF file.
        :type gffs_dir: Path
        """
        dbs_dir = out_dir / "gff_dbs"
        dbs_dir.mkdir(exist_ok=True, parents=True)
        dbfn = dbs_dir / f"{species_prefix}.db"
        gff_file = gffs_dir / f"{species_prefix}.gff"
        if dbfn.exists():
            db = gffutils.FeatureDB(dbfn)
        else:
            db = gffutils.create_db(
                gff_file, dbfn=dbfn, force=True, keep_order=True, merge_strategy="merge"
            )
        with open(gff_file, mode="w") as f:
            for transcript in db.features_of_type("mRNA"):
                raw_id = transcript.attributes.get("ID", [None])[0]
                if raw_id is None:
                    continue
                clean_id = raw_id.split(":")[-1]
                chrom = transcript.chrom.split("_")[-1].replace("Chr", "")
                chrom = str(int(chrom))
                f.write(
                    f"{species_prefix}{chrom}\t{clean_id}\t{transcript.start}\t{transcript.end}\n"
                )

    def prepare_gffs_for_mcscanx(self):
        """
        Prepare all of the original genome GFF files for MCScanX.
        """
        # Fo not need more workers than number of GFF files
        with ProcessPoolExecutor(
            max_workers=min(self.pm.p, len(self.original_gffs))
        ) as executor:
            futures = {
                executor.submit(
                    PangeneConstructor._prepare_individual_gff_for_mcscanx,
                    gff,
                    self.species_to_code_map[gff.parent.name],
                    self.tmp_dir,
                    self.gffs_dir,
                ): gff
                for gff in self.original_gffs
            }
            for future in as_completed(futures):
                gff3 = futures[future]
                self.logger.info(
                    f"Prepared the GFF for {gff3.parent.name} successfully!"
                )

    def _build_blast_dbs(self):
        """
        Creates BLAST databases for every peptide FASTA file.
        """
        blast_dbs_dir = self.tmp_dir / "blast_dbs"
        blast_dbs_dir.mkdir(exist_ok=True, parents=True)
        db_dict = {}
        cmds_dict = {}
        for fasta in self.target_pep_fastas_dir.iterdir():
            genome = fasta.name.replace(".fa", "")
            species_prefix = self.species_to_code_map[genome]
            blast_db_loc = blast_dbs_dir / f"{species_prefix}_prot_db"
            db_dict[species_prefix] = {"db": blast_db_loc, "fasta": fasta}
            cmds_dict[genome] = [
                "makeblastdb",
                "-in",
                fasta,
                "-dbtype",
                "prot",
                "-out",
                blast_db_loc,
            ]

        with ProcessPoolExecutor(
            max_workers=min(self.pm.p, len(self.original_gffs))
        ) as executor:
            futures = {
                executor.submit(
                    execute, cmds, f"Building BLAST database for {genome}"
                ): genome
                for genome, cmds in cmds_dict.items()
            }
            for future in as_completed(futures):
                pass

        self.logger.info(f"Built BLAST databases successfully!")

        return db_dict

    def _run_blastp_queries(self):
        """
        First, builds a BLAST database for each species.
        Then, gets every possible unique pair of species and runs a pairwise query with blastp.
        The parameters used here come in part from the MCScanX documentation suggestions.
        """
        self.logger.info("Building BLAST databases and running pairwise queries.")
        db_dict = self._build_blast_dbs()
        cmds_dict = {}
        for (sp1, d1), (sp2, d2) in product(db_dict.items(), repeat=2):
            db1, fasta1 = d1["db"], d1["fasta"]
            db2, fasta2 = d2["db"], d2["fasta"]
            cmds = [
                "blastp",
                "-db",
                db1,
                "-query",
                fasta2,
                "-evalue",
                "1e-10",
                "-outfmt",
                "6",
                "-max_target_seqs",
                "50",
                "-num_threads",
                "4",
                "-out",
                self.blastps_dir / f"{sp1}_{sp2}.blast",
            ]
            cmds_dict[(sp1, sp2)] = cmds

        with ThreadPoolExecutor(
            max_workers=min(max(1, self.pm.p // 4), len(cmds_dict))
        ) as executor:
            futures = {
                executor.submit(
                    execute, cmds, f"Running blastp for {sp1} querying {sp2}"
                ): (sp1, sp2)
                for (sp1, sp2), cmds in cmds_dict.items()
            }
            for future in as_completed(futures):
                sp1, sp2 = futures[future]

        self.logger.info("Ran pairwise blastp queries successfully!")

    @staticmethod
    def _filter_blast(blast_file: Path, max_diff_r: float):
        """
        Filter out BLAST results such that for each unique query, only results with a score greater than
        or equal to (100 * (1 - `max_diff_r`))% of the top score will be kept.
        `max_diff_r` is 0.001 by default, so kept scores would be 99.9% of the top score for each query.

        :param blast_file: t
        :type blast_file: Path
        :param max_diff_r: The margin to use for query filtering.
        :type max_diff_r: float
        """
        blast_unfiltered = pd.read_csv(
            blast_file,
            sep="\t",
            header=None,
            names=[
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
            ],
        )
        top_scores = blast_unfiltered.groupby("qseqid")["bitscore"].transform("max")
        blast_filtered = blast_unfiltered[
            blast_unfiltered["bitscore"] >= (1.0 - max_diff_r) * top_scores
        ]
        blast_filtered_file = blast_file.with_suffix(".top")
        blast_filtered.to_csv(blast_filtered_file, sep="\t", index=False, header=False)
        # Remove the old file
        blast_file.unlink()

    def blastp_alignment_and_filtering(self):
        """
        Builds blast databases and runs pairwise queries on all species, then filters them to get only top scores.
        """
        self._run_blastp_queries()

        self.logger.info("Filtering BLAST results.")

        blast_files = [b for b in self.blastps_dir.glob("*.blast")]

        self.logger.info(f"Filtering {blast_file.name}")
        with ProcessPoolExecutor(
            max_workers=min(max(1, self.pm.p // 2), len(blast_files))
        ) as executor:
            futures = {
                executor.submit(
                    PangeneConstructor._filter_blast, blast_file, self.max_diff_r
                ): blast_file
                for blast_file in blast_files
            }
            for future in as_completed(futures):
                blast_file = futures[future]

        self.logger.info("Filtered BLAST results successfully!")

    @staticmethod
    def _prep_files_and_run_mcscanx(
        pair: tuple[str, str], mcscanx_dir: Path, blastps_dir: Path, gffs_dir: Path
    ):
        """
        Gathers all BLAST and GFF files into single master files to use for MCScanX, then runs MCScanX.

        :param pair: The species pair to compare.
        :type pair: tuple[str, str]
        :param mcscanx_dir: The directory for MCScanX results and inputs.
        :type mcscanx_dir: Path
        :param blastps_dir: The directory with the BLAST files.
        :type blastps_dir: Path
        :param gffs_dir: The directory with the GFF files.
        :type gffs_dir: Path
        """
        pair_loc = mcscanx_dir / f"{pair[0]}_{pair[1]}"
        blasts = [blastps_dir / f"{'_'.join(p)}.top" for p in product(pair, repeat=2)]
        gffs = [gffs_dir / f"{s}.gff" for s in pair]
        blast_out_file = pair_loc.with_suffix(".blast").resolve()
        gff_out_file = pair_loc.with_suffix(".gff").resolve()
        concat_files(blasts, blast_out_file)
        concat_files(gffs, gff_out_file)

        cmds = ["MCScanX", "-b", "2", pair_loc]
        execute(cmds, "Running MCScanX.")
        blast_out_file.unlink()
        gff_out_file.unlink()

    # TODO: move to utils?
    @staticmethod
    def _split_file(path: Path, chunk_size: int) -> list[Path]:
        """
        Splits a file into chunks of at most `chunk_size`. If the original file is
        file.ext, each split will be named file.part#.ext.

        :param path: The location of the file to split.
        :type path: Path
        :param chunk_size: The desired number of lines to have in each subdivision.
        :type chunk_size: int
        :return: A list of all the split files.
        :rtype: list[Path]
        """
        lines = path.read_text().splitlines(keepends=True)
        chunk_paths = []
        for i in range(0, len(lines), chunk_size):
            chunk_path = path.with_suffix(f".part{i // chunk_size}{path.suffix}")
            chunk_path.write_text("".join(lines[i : i + chunk_size]))
            chunk_paths.append(chunk_path)
        return chunk_paths

    @staticmethod
    def _run_ka_ks_chunk(
        chunk_path: Path, add_k_vals_pl: Path, fasta_reference: Path
    ) -> Path:
        """
        Add Ka and Ks values to the given chunk.

        :param chunk_path: The path to the chunk file.
        :type chunk_path: Path
        :param add_k_vals_pl: The MCScanX script to use to get Ka and Ks values.
        :type add_k_vals_pl: Path
        :param fasta_reference: The FASTA file to reference for calculating Ka and Ks values.
        :type fasta_reference: Path
        :return: The modified chunk file.
        :rtype: Path
        """
        out_path = chunk_path.with_suffix(".col_ks")
        cmds = [
            "perl",
            add_k_vals_pl,
            "-i",
            chunk_path,
            "-d",
            fasta_reference,
            "-o",
            out_path,
        ]
        with tempfile.TemporaryDirectory() as tmp_cwd:
            execute(cmds, f"Finding Ka/Ks for {chunk_path.name}", cwd=tmp_cwd)
        chunk_path.unlink()
        return out_path

    def _add_ka_ks_information(
        self,
        pairs: list[tuple[str, str]],
        add_k_vals_pl: Path,
        fasta_reference: Path,
        min_chunk_size: int = 200,
        oversubscribing_factor: int = 4,
        p: int = 1,
    ):
        """
        Adds Ka and Ks information to each collinearity file output by MCScanX.

        :param pairs: All species pairs.
        :type pairs: list[tuple[str, str]]
        :param add_k_vals_pl: The MCScanX script to use to get Ka and Ks values.
        :type add_k_vals_pl: Path
        :param fasta_reference: The FASTA file to reference for calculating Ka and Ks values.
        :type fasta_reference: Path
        :param min_chunk_size: The minimum size of a chunk to process.
        :type min_chunk_size: int
        :param oversubscribing_factor: How many threads per processor should be used.
        :type oversubscribing_factor: int
        :param p: The number of processors to use.
        :type p: int
        """
        self.logger.info("Adding Ka/Ks information for collinearity files.")
        # Get the .collinearity and .col_ks file locations for each pair
        pair_to_coll_file_dict: dict[tuple[str, str], Path] = {
            pair: (self.mcscanx_dir / f"{pair[0]}_{pair[1]}.collinearity").resolve()
            for pair in pairs
        }
        pair_to_out_file_dict: dict[tuple[str, str], Path] = {
            pair: (self.mcscanx_dir / f"{pair[0]}_{pair[1]}.col_ks").resolve()
            for pair in pairs
        }
        pair_to_chunks_dict: dict[tuple[str, str], list[Path]] = {}

        # Get the total size of all the files to calculate a good chunk size
        total_lines = sum(
            len(coll_file.read_text().splitlines())
            for coll_file in pair_to_coll_file_dict.values()
        )
        chunk_size = max(
            min_chunk_size, np.ceil(total_lines / (p * oversubscribing_factor))
        )

        tasks = []
        # Split all of the collinearity files into chunks to process
        for pair, collinearity_file in pair_to_coll_file_dict.items():
            chunks = PangeneConstructor._split_file(collinearity_file, chunk_size)
            pair_to_chunks_dict[pair] = chunks
            tasks.extend((pair, chunk) for chunk in chunks)

        # Add Ka and Ks values to all chunks
        results: dict[tuple[str, str], list[Path]] = defaultdict(list)
        with ThreadPoolExecutor(max_workers=p) as executor:
            futures = {
                executor.submit(
                    PangeneConstructor._run_ka_ks_chunk,
                    chunk,
                    add_k_vals_pl,
                    fasta_reference,
                ): pair
                for pair, chunk in tasks
            }
            for future in as_completed(futures):
                pair = futures[future]
                results[pair].append(future.result())

        # Recombine all chunks
        for pair, out_chunks in results.items():
            # Sort p0_p1.part#.col_ks by .part# so the output will be in the original order
            out_chunks.sort(
                key=lambda path: int(re.search(r"\.part(\d+)\.", path.name).group(1))
            )
            concat_files(out_chunks, pair_to_out_file_dict[pair])
            for pair in out_chunks:
                pair.unlink()

        self.logger.info(
            "Adding Ka/Ks information for collinearity files successfully!"
        )

    def run_mcscanx_and_get_ks(self):
        """
        Create a master FASTA reference and then run MCScanX for each pair.
        """
        fastas = [
            self.target_cds_fastas_dir / f"{genome}.fa"
            for genome in self.species_to_code_map.keys()
        ]
        combined_cds_file = self.mcscanx_dir / "all_cds.fa"
        concat_files(fastas, combined_cds_file)
        combined_cds_file = combined_cds_file.resolve()

        with ProcessPoolExecutor(
            max_workers=min(max(1, self.pm.p // 4), len(self.pairs))
        ) as executor:
            futures = {
                executor.submit(
                    PangeneConstructor._prep_files_and_run_mcscanx,
                    pair,
                    self.mcscanx_dir,
                    self.blastps_dir,
                    self.gffs_dir,
                ): pair
                for pair in self.pairs
            }
            for future in as_completed(futures):
                self.logger.info(f"MCScanX run for {futures[future]} successfully!")

        self._add_ka_ks_information(
            self.pairs, self.add_k_vals_pl, combined_cds_file, p=self.pm.p
        )

    def _build_synteny_and_tandem_tables(self, conn: sql.Connection):
        """
        Build tables in the pangene SQLite database containing all anchor gene pairs
        from syntenic blocks and all tandem duplicate gene pairs for all species pairs.

        :param conn: The active connection to the SQLite database.
        :type conn: sqlite3.Connection
        """
        for pair in self.code_pairs:
            pair_name = "_".join(pair)
            genomes = [k for k, v in self.species_to_code_map.items() if v in pair]
            genome_pair = frozenset(genomes)
            genome_pair_name = "_".join(sorted(genome_pair))

            col_ks_file = self.mcscanx_dir / f"{pair_name}.col_ks"
            col_ks_df = pd.read_csv(
                col_ks_file,
                sep="\t",
                comment="#",
                header=None,
                names=["block_num", "gene_a", "gene_b", "e_value", "ks", "ka"],
            )
            col_ks_df[["block", "num"]] = col_ks_df["block_num"].str.split(
                "-", expand=True
            )
            col_ks_df["block"] = col_ks_df["block"].str.strip().astype(int)
            col_ks_df["num"] = col_ks_df["num"].str.strip(": ").astype(int)
            col_ks_df = col_ks_df.drop(columns="block_num")
            col_ks_df["genome_pair"] = genome_pair_name
            col_ks_df.to_sql(
                "raw_synteny_blocks", conn, if_exists="append", index=False
            )

            # Build table as well
            tandem_file = self.mcscanx_dir / f"{pair_name}.tandem"
            tandem_pairs = [p.split(",") for p in tandem_file.read_text().splitlines()]
            tandem_df = pd.DataFrame(tandem_pairs, columns=["gene_a", "gene_b"])
            tandem_df.to_sql("tandem_pairs", conn, if_exists="append", index=False)

        conn.execute(f"CREATE INDEX idx_raw_pair ON raw_synteny_blocks(genome_pair);")
        conn.execute(f"CREATE INDEX idx_raw_block ON raw_synteny_blocks(block, ks);")
        conn.commit()

    def _calculate_median_ks_and_filter(self, conn: sql.Connection):
        """
        Find median Ks values for all syntenic anchor pairs. Then, find the mean median Ks value
        (either auto-calculating and finding the maximum mean or using the user-defined genome pairs)
        and filter syntenic blocks out.
        The threshold for filtering syntenic blocks is $\\mu_{K_s} + 3 * \\sigma_{K_s}$.

        :param conn: The active connection to the SQLite database.
        :type conn: sqlite3.Connection
        """
        median_calc_query = """
            CREATE TABLE block_medians AS
                SELECT genome_pair, block, AVG(ks) as median_ks
                FROM (
                    SELECT genome_pair, block, ks,
                        ROW_NUMBER() OVER (PARTITION BY genome_pair, block ORDER BY ks) as rank,
                        COUNT(*) OVER (PARTITION BY genome_pair, block) as total
                    FROM raw_synteny_blocks
                )
                WHERE rank IN ((total + 1) / 2, (total + 2) / 2)
                GROUP BY genome_pair, block;
        """
        conn.execute(median_calc_query)

        medians_df = pd.read_sql_query("SELECT * FROM block_medians", conn)

        ks_threshold_pairs = [
            p for p in ks_threshold_pairs if p in medians_df["genome_pair"].values
        ]
        invalid_threshold_pairs = [
            p for p in ks_threshold_pairs if p not in medians_df["genome_pair"].values
        ]

        for pair in invalid_threshold_pairs:
            self.logger.error(
                f"Reference genome pair for median Ks calculation {pair} is not valid!"
            )

        if len(ks_threshold_pairs) > 0:
            target_medians = medians_df[
                medians_df["genome_pair"].isin(ks_threshold_pairs)
            ]["median_ks"]
            mu = target_medians.mean()
            sigma = target_medians.std()
        else:
            self.logger.info(
                "Calculating maximum mean median Ks value, as no valid reference genome pairs were provided."
            )
            target_pair = (
                medians_df.groupby(by="genome_pair")["median_ks"].mean().idxmax()
            )
            target_medians = medians_df[medians_df["genome_pair"] == target_pair][
                "median_ks"
            ]
            mu = target_medians.mean()
            sigma = target_medians.std()

        ks_threshold = mu + 3 * sigma
        filter_ks_query = f"""
            CREATE TABLE filtered_synteny_blocks AS
                SELECT r.*
                FROM raw_synteny_blocks AS r
                JOIN block_medians AS m 
                    ON r.genome_pair = m.genome_pair AND r.block = m.block
                WHERE m.median_ks < {ks_threshold};
        """
        conn.execute(filter_ks_query)
        self.logger.info("Filtered synteny blocks using Ks values successfully!")

        # Save memory/disk space
        conn.execute("DROP TABLE raw_synteny_blocks;")
        conn.execute("DROP TABLE block_medians;")
        conn.execute("VACUUM;")

        conn.execute(
            "CREATE INDEX idx_filtered_pair ON filtered_synteny_blocks(genome_pair)"
        )
        conn.commit()

    def calculate_ks_threshold_and_filter(self):
        """
        Build the database of microsynteny information and filter it.
        """
        # Build db
        if self.pangene_db.exists():
            self.pangene_db.unlink()
        conn = sql.connect(self.pangene_db)
        # Improve performance (db will likely corrupt if interrupted)
        conn.execute("PRAGMA synchronous = OFF")
        conn.execute("PRAGMA journal_mode = MEMORY")
        conn.commit()

        self._build_synteny_and_tandem_tables(conn)
        self._calculate_median_ks_and_filter(conn)

        conn.close()

    def _subdivide_orthogroup(
        pangene_db: Path, ogid: str, genomes: set[str], row: dict[str, Any]
    ):
        """
        Break an orthologous group into smaller clusters using microsynteny information.

        :param pangene_db: The pangene database to query for microsynteny information.
        :type pangene_db: Path
        :param ogid: The name of the orthologous group.
        :type ogid: str
        :param genomes: The genomes in the orthologous group.
        :type genomes: set[str]
        :param row: The row of the orthologous group file. Each key is a genome (or "OGID") and the values are the genes from that genome in the group (or the OGID).
        :type row: dict[str, Any]
        """
        conn = sql.connect(f"file:{pangene_db}?mode=ro", uri=True, timeout=30.0)
        cursor = conn.cursor()

        gene_genome_pairs = [
            (gene, genome) for genome in genomes for gene in row[genome].split(";")
        ]
        if not gene_genome_pairs:
            return None

        gene_df = pd.DataFrame(gene_genome_pairs, columns=["gene", "genome"])
        genes = list(gene_df["gene"].unique())
        query_params = genes + genes
        placeholder = ",".join(["?"] * len(genes))

        synteny_query = f"""
            SELECT gene_a, gene_b
            FROM filtered_synteny_blocks
            WHERE gene_a IN ({placeholder}) AND gene_b IN ({placeholder});
        """
        tandem_query = f"""
            SELECT gene_a, gene_b
            FROM tandem_pairs
            WHERE gene_a IN ({placeholder}) AND gene_b IN ({placeholder});
        """
        synteny_results = cursor.execute(synteny_query, query_params)
        nodes = synteny_results.fetchall()
        tandem_results = cursor.execute(tandem_query, query_params)
        nodes.extend(tandem_results.fetchall())

        cursor.close()
        conn.close()

        # Final orthogroups are connected components.
        graph: ig.Graph = ig.Graph.TupleList(nodes, directed=False)
        components = graph.connected_components()

        name_to_component_dict = dict(zip(graph.vs["name"], components.membership))
        # If a gene is not present in a connected component, then it will be in the "na" group
        gene_df["OGID"] = gene_df["gene"].apply(
            lambda gene: f"{ogid}_{name_to_component_dict.get(gene, 'na')}"
        )

        new_og_df = gene_df.groupby(by=["OGID", "genome"]).agg(";".join)["gene"]
        new_og_df = new_og_df.unstack().reset_index()
        return new_og_df

    def _stream_og_rows(
        self, og_file: Path
    ) -> Generator[str, set[str], dict[str, Any]]:
        """
        Creates a generator to stream rows from an orthologous group TSV file.

        :param og_file: The orthologous group TSV file.
        :type: Path
        :return: A generator yielding the OGID, the genomes in the orthologous group, and the row data.
        :rtype: Generator[str, set[str], dict[str, Any]]
        """
        with og_file.open(mode="r", encoding="utf-8") as ogf:
            reader = csv.DictReader(ogf, delimiter="\t")
            genomes = set(reader.fieldnames) - {"OGID"}
            for row in reader:
                ogid = row["OGID"]
                yield ogid, genomes, row

    def create_new_orthogroups_with_synteny_info(self, og_of_file: Path) -> Path:
        """
        Takes the orthologous groups created by OrthoFinder and breaks them into smaller
        orthologous groups by clustering on microsynteny information.

        :param og_of_file: The TSV file containing OrthoFinder-created orthologous groups.
        :type og_of_file: Path
        :return: The TSV file containing the smaller orthologous groups.
        :rtype: Path
        """
        # Keep results in a queue to help with streaming
        result_queue = queue.Queue()

        def on_row_complete(future):
            result_queue.put(future)

        new_og_file = self.tmp_dir / "new_ogs.tsv"
        is_first_row_write = True
        row_generator = self._stream_og_rows(og_of_file)
        _, genomes, _ = next(self._stream_og_rows(og_of_file))
        ordered_columns = ["OGID"] + sorted(list(genomes))
        chunk_dfs: list[pd.DataFrame] = []
        buffer_size = self.pm.p * 10
        # It can take a very, very large chunk -- 676 columns, 100 genes per column, 16gb should be about 2300 rows
        many_large_groups_max_chunk_size = 2300
        chunk_size = many_large_groups_max_chunk_size * 32
        active_tasks = 0
        with ProcessPoolExecutor(max_workers=self.pm.p) as executor:
            # Only submit a subset of tasks to save memory
            for _ in range(buffer_size):
                try:
                    ogid, genomes, row = next(row_generator)
                    future = executor.submit(
                        PangeneConstructor._subdivide_orthogroup,
                        self.pangene_db,
                        ogid,
                        genomes,
                        row,
                    )
                    future.add_done_callback(on_row_complete)
                    active_tasks += 1
                except StopIteration:
                    break

            while active_tasks > 0:
                future = result_queue.get()
                active_tasks -= 1
                try:
                    new_og_df = future.result()
                    if new_og_df is not None:
                        chunk_dfs.append(new_og_df)
                except Exception as e:
                    self.logger.error(f"{e}")

                # Add another task when one finishes, if there are any left to add
                try:
                    ogid, genomes, row = next(row_generator)
                    future = executor.submit(
                        PangeneConstructor._subdivide_orthogroup,
                        self.pangene_db,
                        ogid,
                        genomes,
                        row,
                    )
                    future.add_done_callback(on_row_complete)
                    active_tasks += 1
                except StopIteration:
                    pass

                # Write the chunk to the new file
                if len(chunk_dfs) >= chunk_size or (
                    active_tasks == 0 and len(chunk_dfs) > 0
                ):
                    merged_chunk = pd.concat(chunk_dfs, ignore_index=True)
                    merged_chunk = merged_chunk.reindex(columns=ordered_columns)
                    merged_chunk.to_csv(
                        new_og_file,
                        sep="\t",
                        mode="w" if is_first_row_write else "a",
                        header=is_first_row_write,
                        index=False,
                    )
                    is_first_row_write = False
                    chunk_dfs.clear()

        return new_og_file

    def rename_orthogroups_and_remove_empty(self, new_og_file: Path):
        """
        Renames the subdivided orthologous groups and sorts them based on size
        and frequency, removing any empty orthologous groups. Also creates two
        additional tables, `freq_vs_family_size.tsv` and `og_counts.tsv`.

        `freq_vs_family_size.tsv` contains the number of orthologous groups (family size) occurring
        at each frequency, i.e., the number of genomes represented in the group.

        `og_counts.tsv` contains the total number of genomes in the orthologous group (size),
        the number of genomes represented in the group (nonempty_genomes_count), and the classification
        of the group: core, softcore, shell, or private (singleton).

        :param new_og_file: The TSV file containing the smaller orthologous groups.
        :type new_og_file: Path
        """
        new_og_df = pd.read_csv(new_og_file, sep="\t")
        # We don't want to count the ogid
        genome_columns = list(set(new_og_df.columns) - {"OGID"})
        # For each column in the row, the number of genes in that column is 0 if not a string (NaN), or 1 + #_semicolons
        new_og_df["size"] = new_og_df[genome_columns].apply(
            lambda row: sum(
                0 if not isinstance(genes, str) else 1 + genes.count(";")
                for genes in row
            ),
            axis=1,
        )
        new_og_df = new_og_df[new_og_df["size"] > 0]
        new_og_df["nonempty_genomes_count"] = new_og_df[genome_columns].apply(
            lambda row: sum(1 if isinstance(genes, str) else 0 for genes in row), axis=1
        )
        new_og_df = new_og_df.sort_values(
            by=["size", "nonempty_genomes_count", "OGID"],
            ascending=False,
            ignore_index=True,
        )
        max_digits = int(np.ceil(np.log10(new_og_df.index.max())))
        new_og_df["OGID"] = new_og_df.apply(
            lambda row: f"{self.reference}_og{(row.name + 1):0{max_digits}}", axis=1
        )

        num_genomes = len(genome_columns)
        freq_vs_family_size_df = (
            new_og_df.groupby(by=["nonempty_genomes_count"])
            .size()
            .rename("family_size")
            .reset_index()
            .rename(columns={"nonempty_genomes_count": "frequency"})
        )
        freq_vs_family_size_file = self.pangene_dir / "freq_vs_family_size.tsv"
        freq_vs_family_size_df.to_csv(freq_vs_family_size_file, sep="\t", index=False)

        og_counts_df = new_og_df[["OGID", "size", "nonempty_genomes_count"]]
        conditions = [
            og_counts_df["nonempty_genomes_count"] == num_genomes,
            og_counts_df["nonempty_genomes_count"] == (num_genomes - 1),
            (og_counts_df["nonempty_genomes_count"] < (num_genomes - 1))
            & (og_counts_df["nonempty_genomes_count"] > 1),
            og_counts_df["nonempty_genomes_count"] == 1,
        ]
        choices = ["core", "softcore", "shell", "private"]
        og_counts_df["class"] = np.select(conditions, choices, default="unclassified")
        og_counts_file = self.pangene_dir / "og_counts.tsv"
        og_counts_df.to_csv(og_counts_file, sep="\t", index=False)

        sorted_og_df_file = self.pangene_dir / "orthogroups.tsv"
        new_og_df = new_og_df.drop(columns=["size", "nonempty_genomes_count"])
        new_og_df.to_csv(sorted_og_df_file, sep="\t", index=False)

        return sorted_og_df_file

    @staticmethod
    def _extract_gtf(annotation_file: Path, pangene_db: Path, table_name: str):
        """
        Prepares an annotation file (`.gtf` or `.gff`) for use by the pipeline.
        Will gunzip the annotation file and convert it to a `.gtf` file if necessary.

        :param annotation_file: The annotation file to use.
        :type annotation_file: str
        :param pangene_db: The SQLite database to use.
        :type pangene_db: Path
        :param table_name: The name of the table to write to.
        :type table_name: str
        """
        gtfs_temp_dir = tempfile.TemporaryDirectory(suffix="gtf")
        gtfs_dir = Path(gtfs_temp_dir.name)
        annotation_name, annotation_type, is_gzipped = get_name_ext_and_is_gzip(
            annotation_file
        )
        # If not gzipped, this will do nothing
        annotation_file = annotation_file.parent / f"{annotation_name}{annotation_type}"
        if is_gzipped:
            gunzippped_annotation_file = (
                gtfs_dir / f"{annotation_name}{annotation_type}"
            )
            gunzip_file_quiet(annotation_file, gunzippped_annotation_file)
            annotation_file = gunzippped_annotation_file

        converted_to_gtf = False
        # --> GTF if annotation file is GFF
        if "gff" in annotation_type.lower():
            gtf_df = gtfs_dir / f"{annotation_name}.gtf"
            cmds = ["gffread", annotation_file, "-T", "-o", gtf_df]
            execute(cmds)
            converted_to_gtf = True
            annotation_file = gtf_df

        # --> Map file
        if ".gtf" in annotation_type.lower() or converted_to_gtf:
            gtf_df = pd.read_csv(annotation_file, sep="\t", header=None, usecols=[8])
            gtf_df = pd.DataFrame.from_records(
                gtf_df[8]
                .str.split("; ")
                .apply(
                    lambda x: dict(
                        item.strip().split(" ") for item in x if item.strip()
                    )
                )
            )
            gtf_df = gtf_df[["gene_id", "transcript_id"]].replace(
                {'"': "", ";": ""}, regex=True
            )
            gtf_df = gtf_df.drop_duplicates()
            gtf_df = gtf_df.rename(columns={"gene_id": "gene", "transcript_id": "mRNA"})

            conn = sql.connect(pangene_db, timeout=30.0)
            gtf_df.to_sql(table_name, conn, if_exists="append", index=False)
            conn.close()

        else:
            raise NotImplementedError("No valid GFF or GTF file found!")

        gtfs_temp_dir.cleanup()

    def melt_orthogroups(self, sorted_og_df_file: Path) -> Path:
        """
        Creates a melted group file from the sorted and filtered orthologous group file.
        Each row in the melted file represents a single gene and contains:
        * the orthologous group name (OGID)
        * the genome/species the gene belongs to (XGAcc)
        * the transcript id (mRNA)
        * the gene id (gene).

        :param sorted_og_df_file: The sorted and filtered orthologous group TSV file.
        :type sorted_og_df_file: Path
        :return: The location of the melted TSV file.
        :rtype: Path
        """
        gene_mRNA_table_name = "gene_mRNA_table"
        conn = sql.connect(self.pangene_db)
        conn.execute(f"DROP TABLE IF EXISTS {gene_mRNA_table_name}")
        conn.close()
        annotation_files = list(self.pangene_fastas_dir.glob("*/*.g[t|f]f*"))
        gtfs_dir = self.tmp_dir / "gtfs"
        gtfs_dir.mkdir(exist_ok=True, parents=True)
        gene_mRNA_maps: list[pd.DataFrame] = []

        with ProcessPoolExecutor(
            max_workers=min(len(annotation_files), self.pm.p)
        ) as executor:
            futures = [
                executor.submit(
                    PangeneConstructor._extract_gtf,
                    annotation_file,
                    self.pangene_db,
                    gene_mRNA_table_name,
                )
                for annotation_file in annotation_files
            ]
            for future in as_completed(futures):
                pass

        melted_file = self.pangene_dir / "orthogroups.melt"
        og_file = sorted_og_df_file
        original_og_df = pd.read_csv(og_file, sep="\t")

        melted_df = original_og_df.melt(
            id_vars=["OGID"], var_name="XGAcc", value_name="mRNA"
        )
        melted_df["mRNA"] = melted_df["mRNA"].str.split(";")
        melted_df = melted_df.explode(column="mRNA").dropna(ignore_index=True)

        batch_size = 900
        conn = sql.connect(f"file:{self.pangene_db}?mode=ro", uri=True, timeout=30.0)
        cursor = conn.cursor()
        # TODO: Could be parallelized? Maybe not worth it
        for i in range(0, len(melted_df), batch_size):
            mRNA_batch = melted_df.loc[i : (i + batch_size), "mRNA"].values.tolist()
            placeholder = ",".join(["?"] * len(mRNA_batch))
            query = f"""
                SELECT mRNA, gene
                FROM {gene_mRNA_table_name}
                WHERE mRNA IN ({placeholder});
            """
            query_result = cursor.execute(query, mRNA_batch)
            batch_pairs = query_result.fetchall()
            mRNA_to_gene_map = {mRNA: gene for mRNA, gene in batch_pairs}
            melted_df.loc[i : (i + batch_size), "gene"] = melted_df.loc[
                i : (i + batch_size), "mRNA"
            ].apply(lambda mRNA: mRNA_to_gene_map[mRNA])
        cursor.close()
        conn.close()

        melted_df = melted_df.sort_values(by=["OGID", "XGAcc", "mRNA"])
        melted_df.to_csv(melted_file, sep="\t", index=False)

        return melted_file

    def _bgz_genome(self, genome: str):
        """
        Performs a binary gzip of a genome's coding sequence FASTA file.

        :param genome: The name of the genome.
        :type genome: str
        """
        fasta_loc = (
            self.pangene_fastas_dir / genome / f"{genome}_cds.fa.gz"
        )  # TODO: make sure this holds up to snuff -- is this generalized enough?
        bgz_loc = self.tmp_dir / genome / f"{genome}_cds.fa.bgz"
        bgz_loc.parent.mkdir(exist_ok=True, parents=True)
        cmds = ["gzip", "-dc", str(fasta_loc)]
        raw_fasta = subprocess.Popen(cmds, stdout=subprocess.PIPE)
        with open(bgz_loc, mode="wb") as bgz_out:
            cmds = ["bgzip"]
            execute_quiet(cmds, stdin=raw_fasta.stdout, stdout=bgz_out)
        cmds = ["samtools", "faidx", bgz_loc]
        execute_quiet(cmds, f"Building FASTA index file for {genome}.")
        self.logger.log(f"FASTA index for {genome} built successfully!")

    def _extract_data(
        self, genome: str, data: list[tuple[str, str]]
    ) -> dict[str, list[str]]:
        """
        Extracts all coding sequences from the FASTA for a given genome that can be
        found in the the provided data.

        :param genome: The genome whose FASTA will be inspected.
        :type genome: str
        :param data: A list of orthogroup and transcript ids. The transcript ids will be searched for in the FASTA.
        :type data: list[tuple[str, str]]
        :return: All found sequences from the genome. A dictionary mapping all orthogroups to lists of FASTA sequence strings.
        :rtype: dict[str, list[str]]
        """
        fasta_loc = self.tmp_dir / genome / f"{genome}_cds.fa.bgz"
        fasta = Fasta(str(fasta_loc))
        local_buffer = defaultdict(list)
        for orthogroup, mRNA in data:
            if mRNA in fasta:
                sequence = fasta[mRNA][:].seq
                local_buffer[orthogroup].append(f">{mRNA}_{orthogroup}\n{sequence}\n")
            else:
                self.logger.error(f"{mRNA} not found in {genome}")
        return local_buffer

    def prune_redundancies(self):
        """
        Prunes redundancies from the pangene orthogroups using CD-HIT.
        """
        self.logger.info("Shrinking orthogroups by removing redundancies.")
        self.melt_df = pd.read_csv(self.grp_file, sep="\t", dtype=str)
        grouped_melt_df = self.melt_df.groupby("XGAcc")
        extraction_data: dict[str, list[tuple[str, str]]] = {}
        for genome in grouped_melt_df.groups.keys():
            extraction_data[genome] = (
                grouped_melt_df.get_group(genome)[["OGID", "mRNA"]]
                .to_records(index=False)
                .tolist()
            )

        with ThreadPoolExecutor(max_workers=min(self.pm.p + 4, 32)) as executor:
            executor.map(
                self._bgz_genome, [genome for genome in grouped_melt_df.groups.keys()]
            )

        self.logger.info("Extracting FASTA records.")
        global_buffer = defaultdict(list)
        with ProcessPoolExecutor(max_workers=min(16, self.pm.p)) as executor:
            future_to_genome = {
                executor.submit(
                    self._extract_data, genome, extraction_data[genome]
                ): genome
                for genome in extraction_data.keys()
            }
            try:
                for future in tqdm(
                    as_completed(future_to_genome), total=len(future_to_genome)
                ):
                    try:
                        fasta_data = future.result()
                        for orthogroup, sequences in fasta_data.items():
                            global_buffer[orthogroup].extend(sequences)
                    except BrokenProcessPool:
                        self.logger.error(
                            f"A process experienced a critical error while pruning {self.reference}, likely due to too little memory!"
                        )
                        raise
                    except Exception as e:
                        genome = future_to_genome[future]
                        self.logger.error(
                            f"Something went wrong processing {genome} while pruning {self.reference}: {e}"
                        )

            except (KeyboardInterrupt, SystemExit):
                for future in future_to_genome:
                    future.cancel()
                raise

        self.logger.info("FASTA records extracted successfully!")

        fasta_dir = self.pangene_dir / "fastas"
        if fasta_dir.exists():
            shutil.rmtree(fasta_dir)
        fasta_dir.mkdir(exist_ok=True, parents=True)
        full_dir = fasta_dir / "full"
        full_dir.mkdir(exist_ok=True, parents=True)
        self.reduced_dir = fasta_dir / "reduced"
        self.reduced_dir.mkdir(exist_ok=True, parents=True)

        self.logger.info("Writing orthologous groups to FASTA files.")
        for orthogroup, sequences in global_buffer.items():
            out_fasta_loc = (
                full_dir / f"{orthogroup}.fa"
            )  # TODO: consider sharding later
            with open(out_fasta_loc, mode="w") as fout:
                fout.writelines(sequences)
        del global_buffer
        self.logger.info("Orthologous group FASTAs written successfully!")

        self.logger.info("Using CD-HIT to make orthologous groups less redundant.")
        cmds = {
            og: [
                str(cmd)
                for cmd in [
                    "cd-hit-est",
                    "-i",
                    full_dir / f"{og}.fa",
                    "-o",
                    self.reduced_dir / f"{og}.reduced",
                    "-c",
                    self.redunancy_thresh,
                    "-n",
                    self.word_size,
                    "-d",
                    0,
                ]
            ]
            for og in self.melt_df["OGID"].unique()
        }
        with ThreadPoolExecutor(max_workers=min(self.pm.p + 4, 32)) as executor:
            futures = [
                executor.submit(
                    execute_quiet,
                    cmds[og],
                    stdout=subprocess.DEVNULL,
                )
                for og in self.melt_df["OGID"].unique()
            ]
            try:
                for future in tqdm(as_completed(futures), total=len(futures)):
                    pass
            except (KeyboardInterrupt, SystemExit):
                for future in future_to_genome:
                    future.cancel()
                raise
        self.logger.info("Orthologous groups reduced with CD-HIT successfully!")

        self.logger.info("Combining orthologous groups into a single FASTA.")
        cmds = ["find", str(self.reduced_dir), "-name", "*.reduced"]
        reduced = subprocess.Popen(cmds, stdout=subprocess.PIPE)
        cmds = ["xargs", "cat"]
        with gzip.open(self.cds_fasta, mode="wt") as f:
            execute_quiet(
                cmds,
                stdin=reduced.stdout,
                stdout=f,
            )

        self.logger.info("One reduced FASTA file created successfully!")

        annotation_df = self.melt_df[["OGID", "mRNA"]]
        annotation_df["original_transcript_id"] = annotation_df["mRNA"]
        annotation_df["mRNA"] = (
            annotation_df["mRNA"].astype(str) + "_" + annotation_df["OGID"].astype(str)
        )
        annotation_df = annotation_df.rename(
            columns={"OGID": "Geneid", "mRNA": "transcript_id"}
        )
        annotation_df.to_csv(self.annotation_file, sep="\t", index=False)
        self.logger.info("Gene-to-orthologous group map file created successfully!")

    def plot_count_difference(self):
        """
        Creating plots to visualize the results of the redundancy removal.
        """
        self.logger.info(
            "Plotting information about how much redundancy was removed from the pangene."
        )
        original_og_counts = (
            self.melt_df[["OGID", "mRNA"]]
            .groupby(["OGID"])
            .count()
            .rename(columns={"mRNA": "original"})
        )
        clstr_files = list(self.reduced_dir.glob("*.clstr"))
        pruned_og_counts = {}
        for clstr_file in clstr_files:
            og = strip_filename(clstr_file)
            with open(clstr_file, mode="r") as f:
                pruned_og_counts[og] = sum(1 for line in f if ">Cluster" in line)
        pruned_og_counts = pd.DataFrame([pruned_og_counts]).T.rename(
            columns={0: "reduced"}
        )

        counts = pd.merge(
            original_og_counts, pruned_og_counts, left_index=True, right_index=True
        )

        # Original counts histogram
        ax = sns.histplot(data=original_og_counts, x="original", stat="count", bins=100)
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.8)
        ax.set_ylabel("Count")
        ax.set_xlabel("Size of Orthogroup")
        ax.axvline(
            x=original_og_counts["original"].mean(),
            color="r",
            linestyle="--",
            linewidth=1,
            label="Mean Count",
        )
        ax.axvline(
            x=len(self.melt_df["XGAcc"].unique()),
            color="orange",
            linestyle="--",
            linewidth=1,
            label="Number of Genomes",
        )
        ax.set_title("Distribution of Original Orthogroup Sizes")
        ax.legend()
        ax.get_figure().savefig(self.plots_dir / "original_count_hist.png", dpi=600)

        # Reduced counts histogram
        ax = sns.histplot(data=pruned_og_counts, x="reduced", stat="count", bins=100)
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.8)
        ax.set_ylabel("Count")
        ax.set_xlabel("Size of Orthogroup")
        ax.axvline(
            x=pruned_og_counts["reduced"].mean(),
            color="r",
            linestyle="--",
            linewidth=1,
            label="Mean Count",
        )
        ax.axvline(
            x=len(self.melt_df["XGAcc"].unique()),
            color="orange",
            linestyle="--",
            linewidth=1,
            label="Number of Genomes",
        )
        ax.set_title("Distribution of Reduced Orthogroup Sizes")
        ax.legend()
        ax.get_figure().savefig(self.plots_dir / "pruned_count_hist.png", dpi=600)

        # Heatmap for reduction
        x = counts["original"]
        y = counts["reduced"]
        fig, ax = plt.subplots()
        ax.set_xscale("log")
        mesh, x_edges, y_edges, im = ax.hist2d(
            x=x,
            y=y,
            bins=[
                np.logspace(np.log10(x.min()), np.log10(x.max()), 100),
                np.linspace(y.min(), y.max(), 100),
            ],
            cmap="viridis",
            cmin=1,
            norm=colors.LogNorm(vmin=1, vmax=None),
        )
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Number of Groups")
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        line = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), 100)
        ax.plot(
            line,
            line,
            color="r",
            linestyle="--",
            label="No reduction",
            transform=ax.transData,
        )
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel("Original")
        ax.set_ylabel("Reduced")
        ax.legend()
        fig.savefig(self.plots_dir / "reduction_comparison.png", dpi=600)

        self.logger.info("Plots made successfully!")

    def __str__(self) -> str:
        return f"PangeneConstructor {self.reference}"

    def __getstate__(self):
        state = self.__dict__.copy()
        del state["logger"]
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.logger = logging.getLogger(f"{str(self)}")
