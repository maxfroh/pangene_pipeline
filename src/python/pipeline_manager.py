#!/usr/bin/env python3
import shutil
from pathlib import Path
from typing import TypedDict

from .full_pipeline_analysis import FullPipelineAnalyzer
from .pangene_constructor import PangeneConstructor
from .param_manager import ParamManager
from .run_manager import RunManager
from .utils import build_logger


class ConfigDict(TypedDict):
    input: dict[str, str]
    output: dict[str, str]
    dependencies: dict[str, str]
    cores: dict[str, int | bool]
    parameters: dict[str, int | float]
    pangene: dict[str, dict[str, str | float]]
    reference: dict[str, dict[str, str]]
    run: dict[str, dict[str, dict[str, list[str]] | list[str]]]


class PipelineManager:
    def __init__(self, config_dict: ConfigDict):
        self.config_dict = config_dict
        self.output_dir = Path(config_dict["output"]["results_dir"]).resolve()
        # if self.output_dir.exists():
        #     shutil.rmtree(self.output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.pangenes_dir = self.output_dir / "pangenes"
        self.pangenes_dir.mkdir(exist_ok=True, parents=True)
        self.references_dir = self.output_dir / "references"
        self.references_dir.mkdir(exist_ok=True, parents=True)
        self.runs_dir = self.output_dir / "runs"
        self.runs_dir.mkdir(exist_ok=True, parents=True)
        self.logs_dir = self.output_dir / "logs"
        self.logs_dir.mkdir(exist_ok=True, parents=True)
        self.tables_dir = self.output_dir / "tables"
        self.tables_dir.mkdir(exist_ok=True, parents=True)
        self.logger = build_logger(f"{str(self)}")

        params = self.config_dict["cores"] | self.config_dict["parameters"]
        self.pm = ParamManager(**params)

        self.pangenes: dict[str, PangeneConstructor] = {}
        self.references: dict[str, tuple[Path, Path]] = {}
        self.runs: dict[str, RunManager] = {}

        for pangene_name, pangene_info in config_dict["pangene"].items():
            try:
                pangene_info["add_k_vals_pl"] = config_dict["dependencies"][
                    "add_k_vals_pl"
                ]
            except KeyError as ke:
                key = str(ke).strip("'").strip('"')
                self.logger.error(
                    f"Required [dependencies] key {key} is missing from the configuration!"
                )
                raise ke

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
                self.logger.info(f"Making {pc.reference} pangene.")
                pc.construct_pangene()
                self.logger.info(f"{pc.reference} pangene made successfully!")

    def run(self):
        for run_name, run_data in self.config_dict["run"].items():
            run_dir = self.runs_dir / run_name
            run_dir.mkdir(exist_ok=True, parents=True)
            run_pm = self.pm.get_run_variant(**run_data.get("params", {}))
            curr_run = RunManager(
                run_name,
                run_data,
                run_pm,
                self.references_dir,
                run_dir,
                self.pangenes,
                self.references,
            )
            self.runs[run_name] = curr_run

            curr_run.perform_de_analysis()

        fpa = FullPipelineAnalyzer(self.runs, self.pangenes, self.tables_dir)
        fpa.analyze_runs()

    def __str__(self):
        return f"PipelineManager"
