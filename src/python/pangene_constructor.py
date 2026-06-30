#!/usr/bin/env python3
import bisect

from src.python.managers import ManagerDict, ReferenceManager
from src.python.utils import *


class PangeneConstructor:
    def __init__(self, mgr: ManagerDict):
        self.mgr = mgr["manager"]
        self.run = mgr["run"]
        self.reference = mgr["reference"]

    def construct_pangene(self):
        self.get_orthologous_groups_with_orthofinder()
        self.prepare_gffs()
        self.align_cds()
        self.get_syntenic_blocks_with_mcscanx()
        self.get_microsynteny()
        self.separate_orthologous_groups()
        self.assign_representative_genes()
        self.accumulate_final_results()
        self.prune_redundancies()

    def get_orthologous_groups_with_orthofinder(self):
        pass

    def prepare_gffs(self):
        pass

    def align_cds(self):
        pass

    def get_syntenic_blocks_with_mcscanx(self):
        pass

    def get_microsynteny(self):
        pass

    def separate_orthologous_groups(self):
        pass

    def assign_representative_genes(self):
        pass

    def accumulate_final_results(self):
        pass

    def build_pangene_fastas(self):
        pass

    def prune_redundancies(self):
        # get word_size; see https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#user-content-CDHITEST for source of numbers
        thresh = self.mgr.redundancy_thresh
        thresh = max(min(0.75, thresh), 1.0)
        threshholds = [0.8, 0.85, 0.88, 0.9, 0.925, 0.95, 0.975, 1.01]
        word_sizes = [4, 5, 6, 7, 8, 9, 10, 11]
        word_size = word_sizes[bisect.bisect(threshholds, thresh)]

        fasta_file: Path = ""
        reduced_out = fasta_file.parent / (strip_filename(fasta_file) + thresh) + ".fa"

        cmds = [
            "cd-hit-est",
            "-i",
            fasta_file,
            "-o",
            reduced_out,
            "-c",
            thresh,
            "-n",
            word_size,
            "-T",
            self.mgr.p,
        ]

        pass
