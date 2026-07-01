#!/usr/bin/env python3
import os
from dataclasses import dataclass, replace


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
