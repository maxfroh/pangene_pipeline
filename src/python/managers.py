import os

from pathlib import Path
from typing import Union, Any
from types import SimpleNamespace

from .utils import *
from ..default_params import _DEFAULT_PARAMS


class RunManager(SimpleNamespace):
    @staticmethod
    def get_runs_and_references(config: dict[str, Any], default_params: dict[str, Any] = _DEFAULT_PARAMS) -> list[dict[str, Union[str , "RunManager"]]]:
        managers = []
        for run in config["runs"]:
            for reference in config["runs"][run]["references"]:
                managers.append({
                    "run": run,
                    "reference": reference,
                    "manager": RunManager(config, run, reference, default_params)
                })
        return managers
    
    def __init__(self, config: dict[str, Any], run: str, reference: str, defaults: dict):
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
            key = str(ke).strip("\'").strip("\"")
            # add categories to have more information
            logger.error(f"Required [input] key {key} is missing from the configuration!")
            raise ke
        
        try:
            locs["results_dir"] = Path(output_files["results_dir"]) / run
            locs["tmp_dir"] = locs["results_dir"] / "tmp"
            locs["kallisto_dir"] = locs["results_dir"] / "kallisto"
            locs["out_file"] = locs["results_dir"] / output_files.get("out_file", "out.txt")
            locs["err_file"] = locs["results_dir"] / output_files.get("out_file", "err.txt")

        except KeyError as ke:
            key = str(ke).strip("\'").strip("\"")
            # add categories to have more information
            logger.error(f"Required [output] key {key} is missing from the configuration!")
            raise ke
        
        # add additional run information and prune references to other runs
        run_data["run"] = run
        run_data["reference"] = reference
        run_data = {k: run_data[k] for k in ["run", "conditions", "samples"]}
        run_data["sample_names"] = {s: strip_all_extensions(s) for s in run_data["samples"]}
        
        # handle parameters defaulting
        params = config["parameters"]
        for k, v in defaults.items():
            if k not in params:
                params[k] = v
        # determine p!
        if params["auto_allocate_processors"]:
            params["p"] = os.cpu_count() # TODO: 8 or max?????
        
        data = locs | params | run_data
        
        super().__init__(data)
        
        self.create_dirs(data)
        
    def create_dirs(self, data: dict[str, Any]):
        for key in data.keys():
            if "_dir" in key:
                logger.debug(f"Making {data[key]}")
                os.makedirs(data[key], exist_ok=True)
        
    def append_filename(self, original: str | Path, append: str, delimiter: str = "_") -> Path:
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
        if old in self.sample_names:
            self.sample_names[new] = self.sample_names[old]
        self.samples[i] = new
        