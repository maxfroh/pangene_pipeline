#!/usr/bin/env python3
import bisect
import gzip
import logging
import shutil
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pyfaidx import Fasta
from tqdm import tqdm

from .param_manager import ParamManager
from .utils import build_logger, execute, execute_quiet, strip_filename


class PangeneConstructor:
    def __init__(
        self,
        reference: str,
        pangene_dir: Path,
        pangene_info: dict[str, str | float],
        pm: ParamManager,
    ):
        """ """
        self.pm = pm
        self.pangene_dir = pangene_dir / reference
        self.pangene_dir.mkdir(exist_ok=True, parents=True)
        self.tmp_dir = self.pangene_dir / "tmp"
        self.tmp_dir.mkdir(exist_ok=True, parents=True)
        self.plots_dir = self.pangene_dir / "plots"
        self.plots_dir.mkdir(exist_ok=True, parents=True)
        self.reference = reference
        self.constructed = False
        self.logger = build_logger(f"{str(self)}")

        try:
            self.grp_file = Path(pangene_info["grp_file"])
            self.pangene_fastas_dir = Path(pangene_info["pangene_fastas_dir"])
            self.pangene_fastas_dir.mkdir(exist_ok=True, parents=True)
            self.redunancy_thresh = pangene_info["redundancy_thresh"]
        except KeyError as ke:
            key = str(ke).strip("'").strip('"')
            # add categories to have more information
            self.logger.error(
                f"Required [input] key {key} is missing from the configuration!"
            )
            raise ke
        self.annotation_file = self.pangene_dir / "annotation.map"
        self.cds_fasta = self.pangene_dir / f"{reference}_cds.fa.gz"

        if self.redunancy_thresh < 0.75 or self.redunancy_thresh > 1.0:
            thresh = max(min(0.75, self.redunancy_thresh), 1.0)
            self.logger.info(
                f"Pangene redundancy threshold for {self.reference} is out of bounds! New threshold is {thresh} (was {self.redunancy_thresh})"
            )
            thresh = self.redunancy_thresh

        # get word_size; see https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#user-content-CDHITEST for source of numbers
        threshholds = [0.8, 0.85, 0.88, 0.9, 0.925, 0.95, 0.975, 1.01]
        word_sizes = [4, 5, 6, 7, 8, 9, 10, 11]
        self.word_size = word_sizes[bisect.bisect(threshholds, self.redunancy_thresh)]

    def get_reference_info(self) -> tuple[Path, Path]:
        if not self.constructed:
            self.logger.info(f"Pangene {str(self)} not constructed! Constructing now.")
            self.construct_pangene()
        return (self.annotation_file, self.cds_fasta)

    def construct_pangene(self):
        self.get_orthologous_groups_with_orthofinder()
        self.prepare_gffs()
        self.align_cds()
        self.get_syntenic_blocks_with_mcscanx()
        self.get_microsynteny()
        self.separate_orthologous_groups()
        self.assign_representative_genes()
        self.accumulate_final_results()
        # self.prune_redundancies()
        # self.plot_count_difference()
        self.constructed = True

    def get_orthologous_groups_with_orthofinder(self):
        pass

    def prepare_gffs(self):
        pass

    def align_cds(self):
        pass

    def get_syntenic_blocks_with_mcscanx(self):
        pass

    def get_microsynteny(self):
        pass

    def separate_orthologous_groups(self):
        pass

    def assign_representative_genes(self):
        pass

    def accumulate_final_results(self):
        pass

    def build_pangene_fastas(self):
        pass

    def plot_reduction(self):
        pass

    def _bgz_genome(self, genome):
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
            # execute(
            #     cmds,
            #     stdin=raw_fasta.stdout,
            #     stdout=bgz_out,
            #     log_out=False,
            #     log_err=False,
            #     path_out=False,
            #     file_out=True,
            #     logger=self.logger,
            # )
        cmds = ["samtools", "faidx", bgz_loc]
        execute_quiet(cmds, f"Building FASTA index file for {genome}.")
        self.logger.log(f"FASTA index for {genome} built successfully!")

    def _extract_data(self, genome, data: list[tuple[str, str]]):
        fasta_loc = self.tmp_dir / genome / f"{genome}_cds.fa.bgz"
        fasta = Fasta(str(fasta_loc))
        local_buffer = defaultdict(list)
        for orthogroup, mrna in data:
            if mrna in fasta:
                sequence = fasta[mrna][:].seq
                local_buffer[orthogroup].append(f">{mrna}_{orthogroup}\n{sequence}\n")
            else:
                self.logger.error(f"{mrna} not found in {genome}")
        return local_buffer

    def prune_redundancies(self):
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
                    fasta_data = future.result()
                    for orthogroup, sequences in fasta_data.items():
                        global_buffer[orthogroup].extend(sequences)
            except KeyboardInterrupt:
                executor.shutdown(wait=False, cancel_futures=True)
            except MemoryError:  # is this good enough?
                executor.shutdown(wait=False, cancel_futures=True)
                self.logger.error(
                    f"Ran out of memory while extracting FASTA information to build the {self.reference} pangene!"
                )
            finally:
                executor.shutdown(wait=True)
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
            out_fasta_loc = full_dir / f"{orthogroup}.fa"
            with open(out_fasta_loc, mode="w") as fout:
                fout.writelines(sequences)
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
            except KeyboardInterrupt:
                executor.shutdown(wait=False, cancel_futures=True)
            finally:
                executor.shutdown(wait=True)
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

        # for genome in grouped_melt_df.groups.keys():
        shutil.rmtree(self.tmp_dir)

        annotation_df = self.melt_df[["OGID", "mRNA"]]
        annotation_df["mRNA"] = (
            annotation_df["mRNA"].astype(str) + "_" + annotation_df["OGID"].astype(str)
        )
        annotation_df = annotation_df.rename(
            columns={"OGID": "Geneid", "mRNA": "transcript_id"}
        )
        annotation_df.to_csv(self.annotation_file, sep="\t", index=False)
        self.logger.info("Gene-to-orthologous group map file created successfully!")

    def plot_count_difference(self):
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

        # original counts histogram
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

        # reduced counts histogram
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

        # heatmap for reduction
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
