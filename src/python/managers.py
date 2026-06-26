#!/usr/bin/env python3
import os
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Protocol, TypedDict, Union

from src.defaults import _DEFAULT_PARAMS
from src.python.utils import *


class ManagerDict(TypedDict):
    run: str
    reference: str
    manager: "ReferenceManager"


class _Placeholder(SimpleNamespace):
    def __init__(self, config: dict[str, Any]):
        super().__init__(config)


class ReferenceManager(SimpleNamespace):
    @staticmethod
    def get_runs_and_references(
        config: dict[str, Any], default_params: dict[str, Any] = _DEFAULT_PARAMS
    ) -> list[ManagerDict]:
        managers = []
        for run in config["runs"]:
            for reference in config["runs"][run]["references"]:
                manager = ReferenceManager(config, run, reference, default_params)
                manager_dict: ManagerDict = {
                    "run": run,
                    "reference": reference,
                    "manager": manager,
                }
                managers.append(manager_dict)
        return managers

    def __init__(
        self, config: dict[str, Any], run: str, reference: str, defaults: dict
    ):
        self._config = config

        input_files = config["input"]
        output_files = config["output"]
        run_data = config["runs"][run]
        locs = {}

        # build out all files
        try:
            locs["base_dir"] = Path(input_files["base_dir"])
            locs["sample_dir"] = Path(input_files["sample_dir"])
            locs["annotation_file"] = Path(run_data["annotations"][reference])
            locs["reference_file"] = Path(run_data["references"][reference])

        except KeyError as ke:
            key = str(ke).strip("'").strip('"')
            # add categories to have more information
            logger.error(
                f"Required [input] key {key} is missing from the configuration!"
            )
            raise ke

        try:
            locs["results_dir"] = Path(output_files["results_dir"]) / run
            locs["tmp_dir"] = locs["results_dir"] / "tmp"
            locs["kallisto_dir"] = locs["results_dir"] / "kallisto" / reference
            locs["out_file"] = locs["results_dir"] / output_files.get(
                "out_file", "out.txt"
            )
            locs["err_file"] = locs["results_dir"] / output_files.get(
                "out_file", "err.txt"
            )

        except KeyError as ke:
            key = str(ke).strip("'").strip('"')
            # add categories to have more information
            logger.error(
                f"Required [output] key {key} is missing from the configuration!"
            )
            raise ke

        # add additional run information and prune references to other runs
        run_data["run"] = run
        run_data["reference"] = reference
        run_data = {k: run_data[k] for k in ["run", "conditions", "samples"]}
        run_data["sample_name_map"] = {
            s: strip_filename(s) for s in run_data["samples"]
        }

        # handle parameters defaulting
        params = config["parameters"]
        for k, v in defaults.items():
            if k not in params:
                params[k] = v

        # determine threading parameters (total processors, number of subprocesses, actual p per task)
        if params["auto_allocate_processors"]:
            params["total_p"] = os.cpu_count()
        else:
            params["total_p"] = params["p"]
        params["num_threads"] = max(1, params["total_p"] // 8)
        if params["total_p"] < 4:
            params["p"] = 1
        elif 4 <= params["total_p"] < 8:
            params["p"] = 4
        else:
            params["p"] = 8

        data = locs | params | run_data

        super().__init__(data)

        self.create_dirs(data)

    def create_dirs(self, data: dict[str, Any]):
        for key in data.keys():
            if "_dir" in key:
                logger.debug(f"Making {data[key]}")
                os.makedirs(data[key], exist_ok=True)

    def append_filename(
        self, original: str | Path, append: str, delimiter: str = "_"
    ) -> Path:
        full = Path(original)
        # temporarily remove all extensions
        filename, ext = os.path.splitext(full.name)
        filename, ext2 = os.path.splitext(filename)
        if len(ext2) > 0:
            ext = ext2 + ext
        return full.parent / (filename + delimiter + append + ext)

    def update_sample_name(self, i, new):
        # change sample to a new file, but keep its plain name
        old = self.samples[i]
        if old in self.sample_name_map:
            self.sample_name_map[new] = self.sample_name_map[old]
        self.samples[i] = new
