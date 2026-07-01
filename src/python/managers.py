#!/usr/bin/env python3
import os
from copy import deepcopy
from dataclasses import dataclass, field, replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any, NotRequired, Protocol, TypedDict, Union

import pandas as pd

from src.defaults import _DEFAULT_PARAMS
from src.python.deg import DEG
from src.python.pangene_constructor import PangeneConstructor
from src.python.utils import *


class ConfigDict(TypedDict):
    input: dict[str, str]
    output: dict[str, str]
    cores: dict[str, int | bool]
    parameters: dict[str, int | float]
    pangene: dict[str, dict[str, str | float]]
    reference: dict[str, dict[str, str]]
    run: dict[str, dict[str, dict[str, list[str]] | list[str]]]


@dataclass
class ParamManager:
    p = 1  # the total number of cores allocated to the pipeline
    auto_allocate_processors = (
        False  # automatically determine the maximum number of processors available
    )
    alpha = 0.05  # for p-value filtering
    l2FC_thresh = 1  # for l2FC filtering
    frag_length_mean = 200  # for kallisto quant
    frag_length_std = 20  # for kallisto quant
    # redundancy_thresh: 0.98,  # by how much to reduce redunandant genes in the pangene

    def __post_init__(self):
        if self.auto_allocate_processors:
            self.p = os.process_cpu_count()

    def get_run_variant(self, **overrides) -> "ParamManager":
        return replace(self, **overrides)


@dataclass
class ReferenceManager:
    run: str
    reference: str
    pm: ParamManager
    reference_dir: Path
    annotation_file: Path
    cds_fasta: Path
    log_dir: Path
    out_file: Path = field(init=False)
    err_file: Path = field(init=False)
    tmp_dir: Path = field(init=False)
    dge_dir: Path = field(init=False)

    def __post_init__(self):
        self.reference_dir = Path(self.reference_dir)
        self.reference_dir.mkdir(exist_ok=True, parents=True)
        self.cds_fasta = Path(self.cds_fasta)
        self.tmp_dir = self.reference_dir / "tmp"
        self.tmp_dir.mkdir(exist_ok=True, parents=True)
        self.dge_dir = self.reference_dir / "dge"
        self.dge_dir.mkdir(exist_ok=True, parents=True)
        self.out_file = self.log_dir / f"{self.reference}.out"
        self.err_file = self.log_dir / f"{self.reference}.err"
        self.annotation_file = self.prepare_annotation_file(Path(self.annotation_file))

    def prepare_annotation_file(self, annotation_file: Path) -> Path:
        """
        Prepares an annotation file (`.gtf` or `.gff`) for use by the pipeline.
        Will gunzip the annotation file and convert it to a `.gtf` file if necessary.

        :param annotation_file: The annotation file to use.
        :type annotation_file: str
        :return: The location of the usable `.gtf` file.
        :rtype: str
        """
        annotation_name, annotation_type, is_gzipped = get_name_ext_and_is_gzip(
            annotation_file
        )
        annotation_file = annotation_file.parent / f"{annotation_name}{annotation_type}"
        if is_gzipped:
            gunzippped_annotation_file = (
                self.tmp_dir / f"{annotation_name}{annotation_type}"
            )
            gunzip(annotation_file, out_file=gunzippped_annotation_file, replace=False)
            annotation_file = gunzippped_annotation_file

        converted_to_gtf = False
        # --> GTF if annotation file is GFF
        if "gff" in annotation_type.lower():
            gtf = self.tmp_dir / f"{annotation_name}.gtf"
            cmds = ["gffread", annotation_file, "-T", "-o", gtf]
            execute(cmds, f"Converting {annotation_file} to .GTF.")
            default_logger.info("Annotation file prepared successfully!")
            converted_to_gtf = True
            annotation_file = gtf

        # --> Map file
        if ".gtf" in annotation_type.lower() or converted_to_gtf:
            gtf = pd.read_csv(self.annotation_file, sep="\t", header=None, usecols=[8])
            gtf = pd.DataFrame.from_records(
                gtf[8]
                .str.split("; ")
                .apply(
                    lambda x: dict(
                        item.strip().split(" ") for item in x if item.strip()
                    )
                )
            )
            gtf = gtf[["gene_id", "transcript_id"]].replace({'"': "", ";": ""})
            gtf = gtf.rename(columns={"gene_id": "Geneid"})

            annotation_file = self.tmp_dir / f"{annotation_name}.map"
            gtf.to_csv(annotation_file, sep="\t", index=False)

        return Path(annotation_file)

    def __str__(self) -> str:
        return f"ReferenceManager {self.run}::{self.reference}"


@dataclass
class RunManager:
    run: str
    run_data: dict[str, dict[str, list[str]] | list[str]]
    pm: ParamManager
    run_dir: Path
    pangene_dict: dict[str, PangeneConstructor]
    reference_dict: dict[str, tuple[Path, Path]]
    logs_dir: Path = field(init=False)
    tables_dir: Path = field(init=False)
    plots_dir: Path = field(init=False)
    references_dir: Path = field(init=False)
    refms: list[ReferenceManager] = field(default=[], init=False)
    samples: dict[str, Path] = field(default={}, init=False)
    conditions: dict[str, str] = field(default={}, init=False)
    pangene_map_file: Path = field(init=False)

    def __post_init__(self):
        self.logs_dir = self.run_dir / "logs"
        self.logs_dir.mkdir(exist_ok=True, parents=True)
        self.tables_dir = self.run_dir / "tables"
        self.tables_dir.mkdir(exist_ok=True, parents=True)
        self.plots_dir = self.run_dir / "plots"
        self.plots_dir.mkdir(exist_ok=True, parents=True)
        self.references_dir = self.run_dir / "references"
        self.references_dir.mkdir(exist_ok=True, parents=True)

        # get annotation and CDS fasta files for each reference this run
        refm_info: dict[str, tuple[Path, Path]]
        pangene_references: list[str] = self.run_data["use"].get("pangene", [])
        constructed_references: list[str] = self.run_data["use"].get("reference", [])
        if len(pangene_references) > 0:
            for pangene in pangene_references:
                refm_info[pangene] = self.pangene_dict[pangene].get_reference_info()
                self.pangene_map_file = self.pangene_dict[pangene].grp_file
        if len(constructed_references) > 0:
            for reference in constructed_references:
                refm_info[reference] = self.references[reference]

        for reference, (annotation_file, cds_fasta) in refm_info:
            refm = ReferenceManager(
                self.run,
                reference,
                self.pm,
                self.references_dir,
                annotation_file,
                cds_fasta,
                self.logs_dir,
            )
            self.refms.append(refm)

        sample_dir = Path(self.run_data["sample_dir"])
        for sample_file in self.run_data["samples"]:
            sample_path = sample_dir / sample_file
            sample_name = strip_filename(sample_file)
            self.samples[sample_name] = sample_path

    def perform_de_analysis(self):
        deg = DEG(self, self.samples)
        deg.perform_de_analysis(self.refms)

    def __getattr__(self, name):
        if hasattr(self.pm, name):
            return getattr(self.pm, name)
        raise AttributeError(f"{self} does not contain {name}!")

    def __str__(self) -> str:
        return f"RunManager {self.run}"


class PipelineManager:
    def __init__(self, config_dict: ConfigDict):
        self.config_dict = config_dict
        self.output_dir = Path(config_dict["output"]["results_dir"])
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.pangenes_dir = self.output_dir / "pangenes"
        self.pangenes_dir.mkdir(exist_ok=True, parents=True)
        self.runs_dir = self.output_dir / "runs"
        self.runs_dir.mkdir(exist_ok=True, parents=True)
        self.logs_dir = self.output_dir / "logs"
        self.logs_dir.mkdir(exist_ok=True, parents=True)

        params = self.config_dict["cores"] | self.config_dict["parameters"]
        self.pm = ParamManager(**params)
        self.pangenes: dict[str, PangeneConstructor] = {}
        self.references: dict[str, tuple[Path, Path]] = {}
        self.runs: dict[str, RunManager] = {}

        for pangene_name, pangene_info in config_dict["pangene"].items():
            self.pangenes[pangene_name] = PangeneConstructor(
                pangene_name, self.pangenes_dir, pangene_info, self.pm
            )

        for reference_name, reference_info in self.config_dict["reference"].items():
            self.references[reference_name] = (
                Path(reference_info["annotation_file"]),
                Path(reference_info["cds_fasta"]),
            )

    def setup(self):
        for pc in self.pangenes.values():
            if not pc.constructed:
                pc.construct_pangene()

    def run(self):
        for run_name, run_data in self.config_dict["run"].items():
            run_dir = self.runs_dir / run_name
            run_dir.mkdir(exist_ok=True, parents=True)
            run_pm = self.pm.get_run_variant(run_data.get("params", {}))
            curr_run = RunManager(
                run_name, run_data, run_pm, run_dir, self.pangenes, self.references
            )

            curr_run.perform_de_analysis()
