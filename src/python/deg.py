#!/usr/bin/env python3
import shutil
from pathlib import Path

import pandas as pd

from src.python.managers import ReferenceManager, RunManager
from src.python.utils import *


class DEG:
    def __init__(self, runm: RunManager, samples: dict[str, Path]):
        self.runm = runm
        self.samples = samples
        self.column_data_file = self.runm.tables_dir / "column_data.tsv"

    def perform_de_analysis(self, refms: list[ReferenceManager]):
        for refm in refms:
            self.perform_individual_de_analysis(refm)

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
    #     default_logger.info("Reads trimmed successfully!")

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
            refm.kallisto_dir / idx_file,
            "-T",
            refm.tmp_dir,
            "-t",
            self.runm.p,
            refm.cds_fasta,
        ]
        default_logger.debug(f"{idx_file}, {refm.kallisto_dir}")
        default_logger.debug(cmds)

        execute(cmds, "Building kallisto index file.")
        default_logger.info("Kallisto index file built successfully!")
        return idx_file

    def kallisto_quantify(self, refm: ReferenceManager, idx_file: Path):
        """
        Quantifies all samples for a given index file.

        :param idx_file: The kallisto index file to reference.
        :type idx_file: Path
        """
        for item in os.listdir(refm.kallisto_dir):
            if os.path.splitext(item) != ".idx":
                shutil.rmtree(refm.kallisto_dir / item, ignore_errors=True)
        for sample_name, sample_file in self.samples.items():
            cmds = [
                "kallisto",
                "quant",
                "-i",
                refm.kallisto_dir / idx_file,
                "-t",
                self.runm.p,
                "--single",
                "-l",
                self.runm.frag_length_mean,
                "-s",
                self.runm.frag_length_std,
                "-o",
                refm.kallisto_dir / sample_name,
                sample_file,
            ]
            execute(cmds, f"Quantifying {sample_name} with kallisto.")
        default_logger.info("Samples quantified successfully!")

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
        cond_dict = {"sample": [], "condition": []}
        for sample_name in self.samples.keys():
            cond_dict["sample"].append(sample_name)
            cond_dict["condition"].append(
                [cond for cond in conds if cond in sample_name][0]
            )
        cond_df = pd.DataFrame.from_dict(cond_dict)
        cond_df.to_csv(self.column_data_file, sep="\t", index=False)

    def run_deseq(self, refm: ReferenceManager):
        cmds = [
            "Rscript",
            "./src/R/deseq.R",
            refm.kallisto_dir,
            self.column_data_file,
            refm.annotation_file,
            refm.tmp_dir / f"deseq_results_{refm.reference}.tsv",
            self.runm.alpha,
            self.runm.l2FC_thresh,
            *self.samples.keys(),
        ]
        execute(cmds, "Preparing for differential expression analysis with DESeq2.")
        default_logger.info("DESeq2 processing complete!")
