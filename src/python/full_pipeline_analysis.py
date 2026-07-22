#!/usr/bin/env python3
from collections import defaultdict
from functools import reduce
from pathlib import Path

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
from matplotlib_set_diagrams import EulerDiagram

from .pangene_constructor import PangeneConstructor
from .run_manager import RunManager
from .utils import build_logger


class FullPipelineAnalyzer:
    def __init__(self, runs_dict: dict[str, RunManager], pangenes_dict: dict[str, PangeneConstructor], tables_dir: Path):
        self.runs_dict = runs_dict
        self.pangenes_dict = pangenes_dict
        self.tables_dir = tables_dir
        self.reference_tables_dir = self.tables_dir / "references"
        self.reference_tables_dir.mkdir(exist_ok=True, parents=True)

        self.logger = build_logger("Full Pipeline Analyzer")
        
        # all references to keep an eye out for
        self.runs: set[str] = set()
        self.references: set[str] = set()
        for run, rm in runs_dict.items():
            self.runs.add(run)
            self.references = self.references.union({refm.reference for refm in rm.refms})

    def _gather_files(self) -> tuple[dict[str, Path], dict[str, dict[str, list[Path]]]]:
        run_files = {}
        reference_files = {ref: {"counts": [], "abundance": []} for ref in self.references}
        pangene_run_dict = defaultdict(list)
        
        for run, rm in self.runs_dict.items():
            rm = self.runs_dict[run]
            run_files[run] = rm.tables_dir / f"{run}_deg_results.tsv"
            refms = rm.refms
            for refm in refms:
                ref = refm.reference
                dge_dir = refm.dge_dir
                reference_files[ref]["counts"].append(dge_dir / f"counts_{ref}.tsv")
                reference_files[ref]["abundance"].append(dge_dir / f"abundance_{ref}.tsv")
                
            pangene_run_dict[frozenset(rm.pangene_references)].append(run)

        return run_files, reference_files, pangene_run_dict

    def analyze_runs(self):
        """
        Aggregates information from runs into two tables:
        1. A table showing mean TPM and DEG classification for every gene across every reference
        2. A table showing expression fold change and FDR for every gene for every genome comparison across every reference
        
        Aggregates information from references into two tables per reference:
        1. A raw gene count matrix for replicates
        2. A TPM matrix for replicates
        """
        return 
        run_files, reference_files, pangene_run_dict = self._gather_files()
        
        # Aggregate run information
        deg_results_file = self.tables_dir / "deg_results.tsv"
        comparison_table_file = self.tables_dir / "comparisons.tsv"
        
        # group by pangenes used in comparison (likely just one)
        for pangenes, runs in pangene_run_dict.items():
            original_result_dfs: list[pd.DataFrame] = []
            for run in runs:
                run_file = run_files[run]
                deg_results_df = pd.read_csv(run_file, sep="\t")
                columns_to_affix = set(deg_results_df.columns) - {"OGID", "Geneid"}
                deg_results_df = deg_results_df.rename(columns={col: f"{run}_{col}" for col in columns_to_affix})
                str_columns = deg_results_df.select_dtypes(include=["str", "string", "object"]).columns
                original_result_dfs.append(deg_results_df)
        
            master_table = reduce(lambda df1, df2: df1.join(df2, how="outer"), original_result_dfs)
            # Add geneids for missing genes -- pass for now
            # for pangene in pangenes:
            #     map_table = pd.read_csv(self.pangenes_dict[pangene].grp_file, sep="\t", usecols=["OGID", "gene"])
            #     missing_gene_locs = master_table[master_table["Geneid"] == "_MISSING_"].index
            #     master_table.iloc[missing_gene_locs] = np.where(master_table.loc[missing_gene_locs, "OGID"] in map_table["OGID"], map_table["OGID"], "_MISSING_")
            print(master_table)        
        
        # Aggregate reference information
        for ref in self.references:
            count_matrix_file = self.tables_dir / f"{ref}_counts.tsv"
            tpm_matrix_file = self.tables_dir / f"{ref}_tpm.tsv"    
        