#!/usr/bin/env python3
import shutil
from pathlib import Path

import pandas as pd

from src.python.managers import ManagerDict, ReferenceManager
from src.python.utils import *


class DEG:
    def __init__(self, mgr: ManagerDict):
        self.mgr = mgr["manager"]
        self.run = mgr["run"]
        self.reference = mgr["reference"]
        self.column_data_file = self.mgr.tmp_dir / "column_data.tsv"

    def perform_deg(self):
        self.clean_reads()
        self.run_kallisto()
        self.make_condition_table()
        self.perform_deseq_analysis()

    def clean_reads(self, fq_in: Path, fq_out: Path):
        for i in range(len(self.mgr.samples)):
            fq_in = self.mgr.samples[i]
            fq_out = self.mgr.tmp_dir / self.mgr.append_filename(fq_in, "trimmed").name
            cmds = [
                "java",
                "-jar",
                "trimmomatic-0.40.jar",
                "SE",
                "-threads",
                self.mgr.p,
                fq_in,
                fq_out,
                "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10",
                "LEADING:3",
                "TRAILING:3",
                "SLIDINGWINDOW:4:15",
                "MINLEN:36",
            ]
            self.mgr.update_sample_name(i, fq_out)
            execute(cmds, f"Trimming reads for {fq_in}.")
        logger.info("Reads trimmed successfully!")

    def build_kallisto_index(self) -> Path:
        """
        Builds a kallisto index file.

        :return: The location of the index file.
        :rtype: Path
        """

        idx_file = f"{strip_filename(self.mgr.reference_file)}.idx"
        # return idx_file # FLAG
        cmds = [
            "kallisto",
            "index",
            "-i",
            self.mgr.kallisto_dir / idx_file,
            "-T",
            self.mgr.tmp_dir,
            "-t",
            self.mgr.p,
            self.mgr.reference_file,
        ]
        logger.debug(f"{idx_file}, {self.mgr.kallisto_dir}")
        logger.debug(cmds)

        execute(cmds, "Building kallisto index file.")
        logger.info("Kallisto index file built successfully!")
        return idx_file

    def kallisto_quantify(self, idx_file: Path):
        """
        Quantifies all samples for a given index file.

        :param idx_file: The kallisto index file to reference.
        :type idx_file: Path
        """

        for item in os.listdir(self.mgr.kallisto_dir):
            if os.path.splitext(item) != ".idx":
                shutil.rmtree(self.mgr.kallisto_dir / item, ignore_errors=True)
        for sample, sample_name in self.mgr.sample_name_map.items():
            cmds = [
                "kallisto",
                "quant",
                "-i",
                self.mgr.kallisto_dir / idx_file,
                "-t",
                self.mgr.p,
                "--single",
                "-l",
                self.mgr.frag_length_mean,
                "-s",
                self.mgr.frag_length_std,
                "-o",
                self.mgr.kallisto_dir / sample_name,
                self.mgr.sample_dir / sample,
            ]
            execute(cmds, f"Quantifying {sample_name} with kallisto.")
        logger.info("Samples quantified successfully!")

    def run_kallisto(self, idx_file: Path = None):
        """
        Gets abundances for samples with kallisto.

        :param idx_file: The kallisto index file to use. Will create one if not provided.
        :type idx_file: Path
        """

        if idx_file is None:
            idx_file = self.build_kallisto_index()
        self.kallisto_quantify(idx_file)

    def make_condition_table(self):
        conds = self.mgr.conditions
        cond_dict = {"sample": [], "condition": []}
        for sample_name in self.mgr.sample_name_map.values():
            cond_dict["sample"].append(sample_name)
            cond_dict["condition"].append(
                [cond for cond in conds if cond in sample_name][0]
            )
        cond_df = pd.DataFrame.from_dict(cond_dict)
        cond_df.to_csv(self.column_data_file, sep="\t", index=False)

    def perform_deseq_analysis(self):
        cmds = [
            "Rscript",
            "./src/R/deseq.R",
            self.mgr.kallisto_dir,
            self.column_data_file,
            self.mgr.annotation_file,
            self.mgr.tmp_dir / f"deseq_results_{self.reference}.tsv",
            self.mgr.alpha,
            self.mgr.l2FC_thresh,
            *self.mgr.sample_name_map.values(),
        ]
        execute(cmds, "Preparing for differential expression analysis with DESeq2.")
        logger.info("DESeq2 processing complete!")
