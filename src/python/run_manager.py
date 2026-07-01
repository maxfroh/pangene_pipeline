#!/usr/bin/env python3
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from .deg import DEG
from .param_manager import ParamManager
from .reference_manager import ReferenceManager
from .utils import strip_filename

if TYPE_CHECKING:
    from .pangene_constructor import PangeneConstructor
# TODO: LOGGING!


@dataclass
class RunManager:
    run: str
    run_data: dict[str, dict[str, list[str]] | list[str]]
    pm: ParamManager
    run_dir: Path
    pangene_dict: dict[str, PangeneConstructor]
    reference_dict: dict[str, tuple[Path, Path]]
    logger: logging.Logger
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
                refm_info[reference] = self.refms[reference]

        for reference, (annotation_file, cds_fasta) in refm_info:
            refm = ReferenceManager(
                self.run,
                reference,
                self.pm,
                self.references_dir,
                annotation_file,
                cds_fasta,
                self.logs_dir,
                self.logger,
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
