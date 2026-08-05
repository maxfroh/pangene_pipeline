import os
from dataclasses import dataclass, replace


@dataclass
class ParamManager:
    p: int = 1  # the total number of cores allocated to the pipeline
    auto_allocate_processors: bool = (
        False  # automatically determine the maximum number of processors available
    )
    alpha: float = 0.05  # for p-value filtering
    l2FC_thresh: float = 1  # for l2FC filtering
    frag_length_mean: int = 200  # for kallisto quant
    frag_length_std: int = 20  # for kallisto quant
    # redundancy_thresh: 0.98,  # by how much to reduce redunandant genes in the pangene

    def __post_init__(self):
        if self.auto_allocate_processors:
            try:
                self.p = int(len(os.sched_getaffinity(0)) * 0.9)
            except AttributeError:
                self.p = int(os.process_cpu_count() * 0.9)

    def get_run_variant(self, **overrides) -> "ParamManager":
        return replace(self, **overrides)
