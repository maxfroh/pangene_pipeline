#!/usr/bin/env python3
from __future__ import annotations

import os
import shutil
from collections import defaultdict
from concurrent.futures import (ProcessPoolExecutor, ThreadPoolExecutor,
                                as_completed)
from concurrent.futures.process import BrokenProcessPool
from functools import reduce
from itertools import combinations
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
from matplotlib_set_diagrams import EulerDiagram

from .reference_manager import ReferenceManager
from .utils import FixedUpSet, build_logger, execute

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
        self.process_results()

    def perform_individual_de_analysis(self, refm: ReferenceManager):
        # self.clean_reads(refm, "", "") ?
        self.run_kallisto(refm)
        self.make_condition_table()
        self.run_deseq(refm)

    def clean_reads(self, refm: ReferenceManager, fq_in: Path, fq_out: Path):
        for i in range(len(self.samples)):
            cmds = [
                "java",
                "-jar",
                "trimmomatic-0.40.jar",
                "SE",
                "-threads",
                self.runm.p,
                fq_in,
                fq_out,
                "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10",
                "LEADING:3",
                "TRAILING:3",
                "SLIDINGWINDOW:4:15",
                "MINLEN:36",
            ]
            # self.mgr.update_sample_name(i, fq_out)
            execute(cmds, f"Trimming reads for {fq_in}.")
        self.logger.info("Reads trimmed successfully!")

    def kallisto_quantify(self, refm: ReferenceManager, idx_file: Path):
        """
        Quantifies all samples for a given index file.

        :param idx_file: The kallisto index file to reference.
        :type idx_file: Path
        """
        max_p = min(self.runm.p, 8)
        max_workers = min(self.runm.p // max_p, len(self.samples))
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            for sample_name, sample_file in self.samples.items():
                cmds = [
                    "kallisto",
                    "quant",
                    "-i",
                    str(refm.index_file),
                    "-t",
                    max_p,
                    "--single",
                    "-l",
                    self.runm.frag_length_mean,
                    "-s",
                    self.runm.frag_length_std,
                    "-o",
                    str(refm.dge_dir / sample_name),
                    str(sample_file),
                ]
                executor.submit(
                    execute,
                    cmds,
                    f"Quantifying {sample_name} with kallisto using {refm}.",
                )

        self.logger.info("Samples quantified successfully!")

    def run_kallisto(self, refm: ReferenceManager):
        """
        Gets abundances for samples with kallisto.

        :param idx_file: The kallisto index file to use. Will create one if not provided.
        :type idx_file: Path
        """
        self.kallisto_quantify(refm, refm.index_file)

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

    def _get_annotation_file(self, references_dir: Path, ref: str):
        if ref in self.pangene_references:
            map_loc = Path("./results/pangenes") / f"{ref}" / "annotation.map"
        else:
            map_loc = list(Path(references_dir / f"{ref}" / "tmp").glob("*.map"))[0]
        return map_loc

    def process_results(self):
        self.logger.info("Processing results for this run.")
        tables_dir = self.runm.tables_dir
        plots_dir = self.runm.plots_dir
        references_dir = self.runm.references_dir
        refs = [refm.reference for refm in self.runm.refms]
        map_files: dict[str, pd.DataFrame] = {}
        epsilon = 1e-4

        if len(self.pangene_references) > 1:
            raise NotImplementedError

        constructed_refs = list(
            set(refs) - set(self.pangene_references)
        )

        conds_table = pd.read_csv(tables_dir / "column_data.tsv", sep="\t")
        samples = conds_table["sample"].values
        condition_to_samples_map = defaultdict(list)
        for sample in samples:
            condition = conds_table[conds_table["sample"] == sample][
                "condition"
            ].values[0]
            condition_to_samples_map[condition].append(sample)
        condition_pairs = [
            (c1, c2, f"{c2}/{c1}")
            for c1, c2 in combinations(condition_to_samples_map.keys(), 2)
        ]

        for ref in refs:
            specific_ref_dge_dir = references_dir / ref / "dge"
            map_files[ref] = pd.read_csv(
                self._get_annotation_file(references_dir, ref), sep="\t"
            ).drop_duplicates()

        self.logger.info("Building a map file for all references.")
        combined_map: pd.DataFrame = None
        for p_ref in self.pangene_references:
            p_map = map_files[p_ref].set_index("original_transcript_id")
            combined_map = p_map.join(
                [
                    map_files[ref].set_index("transcript_id").add_suffix(f"_{ref}")
                    for ref in constructed_refs
                ],
                how="left",
            )
        combined_map.to_csv(tables_dir / "combined_map.map", sep="\t")
        self.logger.info("Map file built successfully!")

        self.logger.info("Filtering using provided alpha and l2FC threshold.")
        filtered_degs = {}
        for ref in refs:
            # get abundance and padj info
            specific_ref_dge_dir = references_dir / ref / "dge"
            abundance_df = pd.read_csv(
                specific_ref_dge_dir / f"abundance_{ref}.tsv", sep="\t"
            )
            deseq_df = pd.read_csv(
                specific_ref_dge_dir / f"deseq_results_{ref}.tsv", sep="\t"
            )
            abundance_df: pd.DataFrame = abundance_df.set_index("Geneid")
            deseq_df = deseq_df.set_index("Geneid")
            # calculate l2FC for each condition
            for c1, c2, name in condition_pairs:
                c1_avg = abundance_df[condition_to_samples_map[c1]].mean(axis=1)
                c2_avg = abundance_df[condition_to_samples_map[c2]].mean(axis=1)
                abundance_df[name] = np.log2((c2_avg + epsilon) / (c1_avg + epsilon))
            l2fc_names = [name for c1, c2, name in condition_pairs]
            combined_df = abundance_df.join(deseq_df, how="outer")

            # saving just reduced info
            reduced_df = combined_df[["padj", *l2fc_names]]
            reduced_df["noDE"] = reduced_df["padj"] < self.runm.alpha
            reduced_df = reduced_df.drop(columns=["padj"])
            for name in l2fc_names:
                reduced_df[f"{name}_upreg"] = reduced_df[name] > 0
                reduced_df[f"{name}_downreg"] = reduced_df[name] < 0
                reduced_df = reduced_df.drop(columns=[name])

            # filter by p and l2FC thresholds
            filtered_df = combined_df[
                (combined_df["padj"] < self.runm.alpha)
                & (combined_df[l2fc_names].abs() >= self.runm.l2FC_thresh).any(axis=1)
            ]
            for name in l2fc_names:
                filtered_df[f"{name}_upreg"] = filtered_df[name] > 0
                filtered_df[f"{name}_downreg"] = filtered_df[name] > 0
            if ref not in self.pangene_references:
                # add OGID
                curr_ref_map = (
                    combined_map[combined_map[f"Geneid_{ref}"].notna()]
                    .set_index(f"Geneid_{ref}")
                    .rename(columns={"Geneid": "OGID"})["OGID"]
                )
                filtered_df = filtered_df.join(curr_ref_map)
                filtered_df = filtered_df.reset_index().set_index("OGID")
                reduced_df = reduced_df.join(curr_ref_map)
                reduced_df = reduced_df.reset_index().set_index("OGID")
            ref_dir_in_tables_dir = tables_dir / ref
            ref_dir_in_tables_dir.mkdir(exist_ok=True, parents=True)
            reduced_df.to_csv(ref_dir_in_tables_dir / "deg_results.tsv", sep="\t")

            filtered_degs[ref] = filtered_df
        self.logger.info("Results filtered successfully!")

        self.logger.info("Creating an Euler Diagram for this run.")
        fig, ax = plt.subplots(1, 1, layout="constrained")
        EulerDiagram.from_sets(
            sets=[set(df.index) for df in filtered_degs.values()],
            set_labels=filtered_degs.keys(),
            ax=ax,
        )
        fig.suptitle(f"DEGs Detected by Each Reference\n(padj < {self.runm.alpha}; absolute l2FC $\\geq$ {self.runm.l2FC_thresh})")
        fig.savefig(plots_dir / f"{self.runm.run}_venn.png", dpi=600)
        self.logger.info("Euler Diagram created successfully!")
        
        self.logger.info("Results for this run processed successfully!")
