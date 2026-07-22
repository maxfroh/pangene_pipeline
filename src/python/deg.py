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
from upsetplot import from_contents

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
        self.logger = build_logger(f"{self}")

    def perform_de_analysis(self, refms: list[ReferenceManager]):
        for refm in refms:
            self.perform_individual_de_analysis(refm)
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
            map_loc, _ = self.runm.pangene_dict[ref].get_reference_info()
        else:
            map_loc = list(Path(references_dir / f"{ref}" / "tmp").glob("*.map"))[0]
        return map_loc

    def _plot_results(
        self,
        refs: list[str],
        combined_result: pd.DataFrame,
        condition_pairs: tuple[str, str, str],
    ):
        self.logger.info("Creating an Euler Diagram and Upset Plot for this run.")
        plots_dir = self.runm.plots_dir
        euler_results: dict[str, set] = {}
        upset_results: dict[str, set] = {}

        for _, _, pair in condition_pairs:
            for ref in refs:
                euler_results[ref] = set(
                    combined_result[
                        combined_result[f"{ref}_{pair}"].isin(["up", "down"])
                    ]["OGID"]
                )
                up_set = set(
                    combined_result[combined_result[f"{ref}_{pair}"] == "up"]["OGID"]
                )
                down_set = set(
                    combined_result[combined_result[f"{ref}_{pair}"] == "down"]["OGID"]
                )
                upset_results[f"{ref}_up"] = up_set - down_set
                upset_results[f"{ref}_down"] = down_set - up_set
                upset_results[f"{ref}_mixed"] = up_set & down_set

        fig, ax = plt.subplots(1, 1, layout="constrained")
        EulerDiagram.from_sets(
            sets=list(euler_results.values()),
            set_labels=euler_results.keys(),
            ax=ax,
        )
        fig.suptitle(
            f"DEGs Captured by Each Reference\n(padj < {self.runm.alpha}; absolute l2FC $\\geq$ {self.runm.l2FC_thresh})"
        )
        fig.savefig(plots_dir / f"{self.runm.run}_venn.png", dpi=600)
        plt.cla()
        self.logger.info("Euler Diagram created successfully!")

        upset_data = from_contents(upset_results)
        category_order = sorted(
            upset_data.index.names,
            key=lambda s: (
                s.rsplit("_")[-1],
                s.removesuffix("_up").removesuffix("_down").removesuffix("_mixed"),
            ),
        )
        upset_data = upset_data.reorder_levels(category_order)
        threshold = 0.01
        upset_data = (
            upset_data.groupby(level=upset_data.index.names)
            .size()
            .astype(int)
            .sort_values(ascending=False)
        )
        large_categories = list(
            upset_data[upset_data / upset_data.sum() > threshold].index
        )
        upset_data = upset_data[upset_data.index.isin(large_categories)]
        levels_to_drop = [
            index
            for index in upset_data.index.names
            if set(upset_data.index.get_level_values(index)) == {False}
        ]
        upset_data = upset_data.droplevel(levels_to_drop)

        upset = FixedUpSet(
            data=upset_data,
            sort_by="-cardinality",
            sort_categories_by="input",
            min_subset_size="1%",
            show_counts=True,
            element_size=40,
        )
        axes = upset.plot()
        intersections_ax = axes["intersections"]
        current_xmin, current_xmax = intersections_ax.get_xlim()
        intersections_ax.margins(y=0.2)
        intersections_ax.set_xlim(current_xmin, current_xmax * 1.01)
        plt.suptitle(
            f"Intersections of OGID up/downregulation\n(Only subsets with size $\\geq$ {threshold} shown)"
        )
        plt.savefig(plots_dir / f"{self.runm.run}_upset.png", dpi=600)
        plt.cla()
        self.logger.info("Upset Plot created successfully!")

    def process_results(self):
        self.logger.info("Processing results for this run.")
        tables_dir = self.runm.tables_dir
        references_dir = self.runm.references_dir
        refs = [refm.reference for refm in self.runm.refms]
        map_files: dict[str, pd.DataFrame] = {}
        epsilon = 1e-4

        if len(self.pangene_references) > 1:
            raise NotImplementedError

        constructed_refs = list(set(refs) - set(self.pangene_references))

        conds_table = pd.read_csv(tables_dir / "column_data.tsv", sep="\t")
        samples = conds_table["sample"].values
        condition_to_samples_map = defaultdict(list)
        for sample in samples:
            condition = conds_table[conds_table["sample"] == sample][
                "condition"
            ].values[0]
            condition_to_samples_map[condition].append(sample)
        condition_pairs = [
            (c1, c2, f"{c1}_{c2}")
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
        combined_results: dict[str, pd.DataFrame] = {}

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
                abundance_df[f"{c1}_mean_TPM"] = c1_avg
                c2_avg = abundance_df[condition_to_samples_map[c2]].mean(axis=1)
                abundance_df[f"{c2}_mean_TPM"] = c2_avg
                abundance_df[f"{name}_FC"] = np.where(
                    c1_avg != 0, c2_avg / c1_avg.replace(0, 1), 0
                )
                abundance_df[f"{name}_l2FC"] = np.where(
                    (c1_avg != 0) | (c2_avg != 0), np.log2((c2_avg / c1_avg.replace(0, 1)).replace(0, 1)), 0
                )
            l2FC_names = [name for c1, c2, name in condition_pairs]
            pair_columns = [
                col
                for c1, c2, name in condition_pairs
                for col in (
                    name,
                    f"{name}_l2FC",
                    f"{name}_FC",
                    f"{c1}_mean_TPM",
                    f"{c2}_mean_TPM",
                )
            ]
            combined_df = abundance_df.join(deseq_df, how="outer")

            # saving just reduced info
            for name in l2FC_names:
                is_de = (combined_df["padj"] < self.runm.alpha) & (
                    combined_df[f"{name}_l2FC"].abs() >= self.runm.l2FC_thresh
                )
                conditions = [
                    ~is_de,
                    (is_de & (combined_df[f"{name}_l2FC"] > 0)),
                    (is_de & (combined_df[f"{name}_l2FC"] < 0)),
                ]
                choices = ["no", "up", "down"]
                combined_df[name] = np.select(conditions, choices, "no")
            combined_df = combined_df[["padj", *pair_columns]]
            if ref not in self.pangene_references:
                # add OGID
                curr_ref_map = (
                    combined_map[combined_map[f"Geneid_{ref}"].notna()]
                    .set_index(f"Geneid_{ref}")
                    .rename(columns={"Geneid": "OGID"})["OGID"]
                )
                combined_df = combined_df.join(curr_ref_map)
                combined_df = combined_df.reset_index()
                # combined_df = combined_df.set_index("OGID", append=True).swaplevel(0, 1)
            else:
                combined_df = combined_df.rename_axis("OGID").reset_index()
            combined_df = combined_df.rename(
                columns={
                    c: f"{ref}_{c}"
                    for c in combined_df.columns
                    if c not in ["Geneid", "OGID"]
                }
            )
            combined_results[ref] = combined_df

            def merge_on_indices_dynamic(left, right):
                if "Geneid" in left.columns and "Geneid" in right.columns:
                    keys = ["OGID", "Geneid"]
                else:
                    keys = ["OGID"]
                return pd.merge(left, right, on=keys, how="outer")

        combined_result = reduce(merge_on_indices_dynamic, combined_results.values())
        target_cols = ["OGID", "Geneid"]
        combined_result = combined_result[
            target_cols
            + [col for col in combined_result.columns if col not in target_cols]
        ]

        self.logger.info("Results filtered successfully!")

        self._plot_results(refs, combined_result, condition_pairs)

        str_columns = ["Geneid"] + [f"{ref}_{name}" for c1, c2, name in condition_pairs for ref in refs]
        print(str_columns)
        # turn NA strings into None to prevent read warnings later on
        combined_result[str_columns] = combined_result[str_columns].fillna("_MISSING_")
        combined_result.to_csv(
            tables_dir / f"{self.runm.run}_deg_results.tsv", sep="\t", index=False
        )

        self.logger.info("Results for this run processed successfully!")

    def __str__(self):
        return f"DEG Analyzer for {self.runm}"
