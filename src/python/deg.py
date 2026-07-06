#!/usr/bin/env python3
from __future__ import annotations

import os
import shutil
from functools import reduce
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib_set_diagrams import EulerDiagram
from upsetplot import UpSet
from itertools import combinations
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns

from .reference_manager import ReferenceManager
from .utils import build_logger, execute

if TYPE_CHECKING:
    from run_manager import RunManager


class DEG:
    def __init__(
        self, runm: RunManager, samples: dict[str, Path], pangene_references: list[str]
    ):
        self.runm = runm
        self.samples = samples
        self.pangene_references = pangene_references
        self.column_data_file = self.runm.tables_dir / "column_data.tsv"
        self.logger = build_logger(f"DEG for {self.runm}")

    def perform_de_analysis(self, refms: list[ReferenceManager]):
        for refm in refms:
            # self.perform_individual_de_analysis(refm)
            pass
        self.build_tables()
        self.plot_results()

    def perform_individual_de_analysis(self, refm: ReferenceManager):
        # self.clean_reads()
        # self.run_kallisto(refm)
        self.make_condition_table()
        self.run_deseq(refm)

    # def clean_reads(self, fq_in: Path, fq_out: Path):
    #     for i in range(len(self.mgr.samples)):
    #         fq_in = self.mgr.samples[i]
    #         fq_out = self.mgr.tmp_dir / self.mgr.append_filename(fq_in, "trimmed").name
    #         cmds = [
    #             "java",
    #             "-jar",
    #             "trimmomatic-0.40.jar",
    #             "SE",
    #             "-threads",
    #             self.mgr.p,
    #             fq_in,
    #             fq_out,
    #             "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10",
    #             "LEADING:3",
    #             "TRAILING:3",
    #             "SLIDINGWINDOW:4:15",
    #             "MINLEN:36",
    #         ]
    #         self.mgr.update_sample_name(i, fq_out)
    #         execute(cmds, f"Trimming reads for {fq_in}.")
    #     self.logger.info("Reads trimmed successfully!")

    def build_kallisto_index(self, refm: ReferenceManager) -> Path:
        """
        Builds a kallisto index file.

        :return: The location of the index file.
        :rtype: Path
        """
        idx_file = f"{refm.reference}.idx"
        # return idx_file # FLAG
        cmds = [
            "kallisto",
            "index",
            "-i",
            refm.dge_dir / idx_file,
            "-T",
            refm.tmp_dir,
            "-t",
            min(self.runm.p, 8),
            refm.cds_fasta,
        ]
        self.logger.debug(f"{idx_file}, {refm.dge_dir}")
        self.logger.debug(cmds)

        execute(cmds, f"Building kallisto index file for {refm}.")
        self.logger.info("Kallisto index file built successfully!")
        return idx_file

    def kallisto_quantify(self, refm: ReferenceManager, idx_file: Path):
        """
        Quantifies all samples for a given index file.

        :param idx_file: The kallisto index file to reference.
        :type idx_file: Path
        """
        for item in os.listdir(refm.dge_dir):
            if os.path.splitext(item) != ".idx":
                shutil.rmtree(refm.dge_dir / item, ignore_errors=True)
        for sample_name, sample_file in self.samples.items():
            cmds = [
                "kallisto",
                "quant",
                "-i",
                refm.dge_dir / idx_file,
                "-t",
                min(self.runm.p, 8),
                "--single",
                "-l",
                self.runm.frag_length_mean,
                "-s",
                self.runm.frag_length_std,
                "-o",
                refm.dge_dir / sample_name,
                sample_file,
            ]
            execute(cmds, f"Quantifying {sample_name} with kallisto using {refm}.")
        self.logger.info("Samples quantified successfully!")

    def run_kallisto(self, refm: ReferenceManager, idx_file: Path = None):
        """
        Gets abundances for samples with kallisto.

        :param idx_file: The kallisto index file to use. Will create one if not provided.
        :type idx_file: Path
        """
        if idx_file is None:
            idx_file = self.build_kallisto_index(refm)
        self.kallisto_quantify(refm, idx_file)

    def make_condition_table(self):
        conds = self.runm.conditions
        print(conds)
        cond_dict = {"sample": [], "condition": []}
        for sample_name in self.samples.keys():
            cond_dict["sample"].append(sample_name)
            cond_dict["condition"].append(
                [cond for cond in conds if cond in sample_name][0]
            )
        cond_df = pd.DataFrame.from_dict(cond_dict)
        cond_df.to_csv(self.column_data_file, sep="\t", index=False)

    def _get_deseq_results_file(self, refm: ReferenceManager) -> str:
        return refm.dge_dir / f"deseq_results_{refm.reference}.tsv"

    def _get_kallisto_counts_file(self, refm: ReferenceManager) -> str:
        return refm.dge_dir / f"counts_{refm.reference}.tsv"

    def _get_kallisto_abundance_file(self, refm: ReferenceManager) -> str:
        return refm.dge_dir / f"abundance_{refm.reference}.tsv"

    def run_deseq(self, refm: ReferenceManager):
        cmds = [
            "Rscript",
            "./src/R/deseq.R",
            refm.dge_dir,
            self.column_data_file,
            refm.annotation_file,
            self._get_deseq_results_file(refm),
            self._get_kallisto_counts_file(refm),
            self._get_kallisto_abundance_file(refm),
            *self.samples.keys(),
        ]
        execute(cmds, "Preparing for differential expression analysis with DESeq2.")
        self.logger.info("DESeq2 processing complete!")

    def _get_annotation_file(self, references_dir: Path, refm: str):
        if refm in self.pangene_references:
            map_loc = Path("./results/pangenes") / f"{refm}" / "annotation.map"
        else:
            map_loc = list(Path(references_dir / f"{refm}" / "tmp").glob("*.map"))[0]
        return map_loc

    def build_tables(self):
        tables_dir = self.runm.tables_dir
        references_dir = self.runm.references_dir
        abundance_results: dict[str, pd.DataFrame] = {}
        counts_results: dict[str, pd.DataFrame] = {}
        deseq_results: dict[str, pd.DataFrame] = {}
        map_files: dict[str, pd.DataFrame] = {}
        combined_map: pd.DataFrame = None

        self.logger.info("Reading results from kallisto and DESeq2.")
        for refm in self.runm.refms:
            ref = refm.reference
            specific_ref_dge_dir = references_dir / ref / "dge"

            map_files[ref] = pd.read_csv(
                self._get_annotation_file(references_dir, ref), sep="\t"
            ).drop_duplicates()

            abundance_results[ref] = pd.read_csv(
                specific_ref_dge_dir / f"abundance_{ref}.tsv", sep="\t"
            )
            counts_results[ref] = pd.read_csv(
                specific_ref_dge_dir / f"counts_{ref}.tsv", sep="\t"
            )
            deseq_results[ref] = pd.read_csv(
                specific_ref_dge_dir / f"deseq_results_{ref}.tsv", sep="\t"
            )

            # if not a pangene, convert genes to transcript
            if "original_transcript_id" in map_files[ref].columns:
                map_files[ref] = map_files[ref].rename(
                    columns={
                        "original_transcript_id": "transcript_id",
                        "transcript_id": "original_transcript_id",
                        "Geneid": "OGID",
                    }
                )

            map_files[ref] = map_files[ref].set_index("transcript_id")

        # build main map
        pangene_maps = {}
        for p_ref in self.pangene_references:
            p_map = map_files[p_ref]
            combined_map = p_map.join(
                [
                    map_files[ref].add_suffix(f"_{ref}")
                    for ref in set([refm.reference for refm in self.runm.refms])
                    - {p_ref}
                ],
                how="left",
            )

        # Only support single pangene comparisons right now
        # if len(pangene_refs) > 1:
        #     combined_map = pangene_maps[pangene_refs[0]]
        #     for p_ref in pangene_refs[1:]:
        #         combined_map = combined_map.join(pangene_maps[p_ref], lsuffix=f"_{pangene_refs[0]}", rsuffix=f"_{p_ref}")
        # else:
        #     combined_map = list(pangene_maps.values())[0]
        combined_map.to_csv(tables_dir / "combined_map.map", sep="\t")

        self.logger.info("Annotating results for easier comparison.")
        for refm in self.runm.refms:
            ref = refm.reference
            self.logger.debug(f"{ref}")
            if ref not in self.pangene_references:
                abundance_results[ref] = abundance_results[ref].merge(
                    map_files[ref].reset_index(), on="Geneid", how="left"
                )
                abundance_results[ref]["OGID"] = combined_map.loc[
                    abundance_results[ref]["transcript_id"]
                ]["OGID"].values
                abundance_results[ref] = abundance_results[ref].set_index("OGID")
                counts_results[ref] = counts_results[ref].merge(
                    map_files[ref].reset_index(), on="Geneid", how="left"
                )
                counts_results[ref]["OGID"] = combined_map.loc[
                    counts_results[ref]["transcript_id"]
                ]["OGID"].values
                counts_results[ref] = counts_results[ref].set_index("OGID")
                deseq_results[ref] = deseq_results[ref].merge(
                    map_files[ref].reset_index(), on="Geneid", how="left"
                )
                deseq_results[ref]["OGID"] = combined_map.loc[
                    deseq_results[ref]["transcript_id"]
                ]["OGID"].values
                deseq_results[ref] = deseq_results[ref].set_index("OGID")
            else:
                abundance_results[ref] = abundance_results[ref].rename(
                    columns={"Geneid": "OGID"}
                )
                counts_results[ref] = counts_results[ref].rename(
                    columns={"Geneid": "OGID"}
                )
                deseq_results[ref] = deseq_results[ref].rename(
                    columns={"Geneid": "OGID"}
                )
                abundance_results[ref] = abundance_results[ref].set_index("OGID")
                counts_results[ref] = counts_results[ref].set_index("OGID")
                deseq_results[ref] = deseq_results[ref].set_index("OGID")

        for refm in self.runm.refms:
            ref = refm.reference
            ref_table_dir = tables_dir / ref
            ref_table_dir.mkdir(parents=True, exist_ok=True)
            abundance_results[ref].to_csv(ref_table_dir / "abundance.tsv", sep="\t")
            counts_results[ref].to_csv(ref_table_dir / "counts.tsv", sep="\t")
            deseq_results[ref].to_csv(ref_table_dir / "deseq.tsv", sep="\t")

        self.logger.info("Annotated tables created and saved successfully!")

    def plot_results(self):
        tables_dir = self.runm.tables_dir
        plots_dir = self.runm.plots_dir

        abundance_results = {}
        counts_results = {}
        deseq_results = {}
        for ref, subdir, file in tables_dir.walk():
            if ref != tables_dir:
                abundance_results[ref.name] = pd.read_csv(
                    ref / "abundance.tsv", sep="\t"
                )
                counts_results[ref.name] = pd.read_csv(
                    ref / "counts.tsv", sep="\t"
                )
                deseq_results[ref.name] = pd.read_csv(
                    ref / "deseq.tsv", sep="\t"
                )

        # Euler Diagram for just the raw counts
        shared_ogs = {
            ref: set(counts_results[ref]["OGID"].values) for ref in counts_results
        }
        fig, ax = plt.subplots(1, 1, layout="constrained")
        EulerDiagram.from_sets(
            sets=list(shared_ogs.values()),
            set_labels=list(shared_ogs.keys()),
            ax=ax,
            cost_function_objective="logarithmic",
        )
        fig.suptitle("Counts of Expressed Genes", fontweight="bold")
        fig.subplots_adjust(top=0.925)
        fig.savefig(plots_dir / "raw_counts_venn.png", dpi=600)
        
        # PAV with TPM
        deg_results = {}
        epsilon = 1e-4
        conds_table = pd.read_csv(tables_dir / "column_data.tsv", sep="\t")
        samples = conds_table["sample"].values
        conditions_to_sample_map = defaultdict(list)
        for sample in samples:
            condition = conds_table[conds_table["sample"] == sample]["condition"].values[0]
            conditions_to_sample_map[condition].append(sample)
            pairs = []
        condition_pairs = [(c1, c2, f"{c1}:{c2}") for c1, c2 in combinations(conditions_to_sample_map.keys(), 2)]
        ref_pairs = [(r1, r2, f"{r1} vs. {r2}") for r1, r2 in combinations(deseq_results.keys(), 2)]

        for ref in deseq_results:
            deg = deseq_results[ref][deseq_results[ref]["padj"] < self.runm.alpha]
            if ref in ["Pangene"]:
                deg = deg.set_index("OGID")
                deg = deg.join(abundance_results[ref].set_index("OGID"))
            else:
                deg = deg.set_index("transcript_id")[["padj"]]
                deg = deg.join(abundance_results[ref].set_index("transcript_id"))
            deg = deg.reset_index()
            
            for condition in conditions_to_sample_map.keys():
                deg[condition] = deg[conditions_to_sample_map[condition]].mean(axis=1)
            
            # l2FC
            for c1, c2, name in condition_pairs:
                deg[name] = np.log2((deg[c1] + epsilon) / (deg[c2] + epsilon))
            for c1, c2, name in condition_pairs:
                deg = deg[np.abs(deg[name]) > self.runm.l2fc_thresh]
                deg[f"{name}_upreg"] = deg[name] > 0
                
            deg_results[ref] = deg
            deg.to_csv(tables_dir / ref / "deg.tsv", sep="\t")
            
        
            