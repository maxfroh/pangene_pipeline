#!/usr/bin/env python3
from __future__ import annotations

import os
import shutil
from functools import reduce
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from .reference_manager import ReferenceManager
from .utils import build_logger, execute

if TYPE_CHECKING:
    from run_manager import RunManager


class DEG:
    def __init__(self, runm: RunManager, samples: dict[str, Path]):
        self.runm = runm
        self.samples = samples
        self.column_data_file = self.runm.tables_dir / "column_data.tsv"
        self.logger = build_logger(f"DEG for {self.runm}")

    def perform_de_analysis(self, refms: list[ReferenceManager]):
        print(refms)
        for refm in refms:
            print(f"Looking at {refm}")
            self.perform_individual_de_analysis(refm)
        print("Done")
        self.build_tables()
        self.plot_results()

    def perform_individual_de_analysis(self, refm: ReferenceManager):
        # self.clean_reads()
        self.run_kallisto(refm)
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

    def build_tables(self):
        tables_dir = self.runm.tables_dir
        references_dir = self.runm.references_dir

        deseq_results: dict[str, pd.DataFrame] = {}
        counts_results: dict[str, pd.DataFrame] = {}
        abundance_results: dict[str, pd.DataFrame] = {}
        map_files: dict[str, pd.DataFrame] = {}
        for refm in self.runm.refms:
            results_file = self._get_deseq_results_file(refm)
            deseq_result = pd.read_csv(results_file, sep="\t", index_col=0)
            counts_file = self._get_kallisto_counts_file(refm)
            counts_result = pd.read_csv(counts_file, sep="\t", index_col=0)
            abundance_file = self._get_kallisto_abundance_file(refm)
            abundance_result = pd.read_csv(abundance_file, sep="\t", index_col=0)
            map_df = pd.read_csv(refm.annotation_file, sep="\t")

    def plot_results(self):
        tables_dir = self.runm.tables_dir
        plots_dir = self.runm.plots_dir
