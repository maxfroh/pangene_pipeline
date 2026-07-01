#!/usr/bin/env python3
from pathlib import Path
from typing import TypedDict

from .pangene_constructor import PangeneConstructor
from .param_manager import ParamManager
from .run_manager import RunManager


class ConfigDict(TypedDict):
    input: dict[str, str]
    output: dict[str, str]
    cores: dict[str, int | bool]
    parameters: dict[str, int | float]
    pangene: dict[str, dict[str, str | float]]
    reference: dict[str, dict[str, str]]
    run: dict[str, dict[str, dict[str, list[str]] | list[str]]]


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
                run_name,
                run_data,
                run_pm,
                run_dir,
                self.pangenes,
                self.references,
            )

            curr_run.perform_de_analysis()
