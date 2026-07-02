#!/usr/bin/env python3
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from .utils import build_logger, execute, get_name_ext_and_is_gzip, gunzip

if TYPE_CHECKING:
    from .param_manager import ParamManager


@dataclass
class ReferenceManager:
    run: str
    reference: str
    pm: ParamManager
    reference_dir: Path
    annotation_file: Path
    cds_fasta: Path
    log_dir: Path
    logger: logging.Logger = field(init=False)
    out_file: Path = field(init=False)
    err_file: Path = field(init=False)
    tmp_dir: Path = field(init=False)
    dge_dir: Path = field(init=False)

    def __post_init__(self):
        self.logger = build_logger(f"{str(self)}")
        self.reference_dir = self.reference_dir / self.reference
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
            os.remove(annotation_file)
            self.logger.info("Annotation file prepared successfully!")
            converted_to_gtf = True
            annotation_file = gtf

        # --> Map file
        if ".gtf" in annotation_type.lower() or converted_to_gtf:
            gtf = pd.read_csv(annotation_file, sep="\t", header=None, usecols=[8])
            gtf = pd.DataFrame.from_records(
                gtf[8]
                .str.split("; ")
                .apply(
                    lambda x: dict(
                        item.strip().split(" ") for item in x if item.strip()
                    )
                )
            )
            gtf = gtf[["gene_id", "transcript_id"]].replace(
                {'"': "", ";": ""}, regex=True
            )
            gtf = gtf.rename(columns={"gene_id": "Geneid"})

            annotation_file = self.tmp_dir / f"{annotation_name}.map"
            gtf.to_csv(annotation_file, sep="\t", index=False)

        return Path(annotation_file)

    def __str__(self) -> str:
        return f"ReferenceManager {self.run}::{self.reference}"
