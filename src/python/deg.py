import shutil

from pathlib import Path

from .managers import RunManager
from .utils import *


class DEG():
    def __init__(self, mgr: RunManager):
        self.mgr = mgr
        
    def perform_deg(self):
        # self.clean_reads()
        self.run_kallisto()
        self.scale_abundances()
    
    def clean_reads(self, fq_in: Path, fq_out: Path):
        for i in range(len(self.mgr.samples)):
            fq_in = self.mgr.samples[i]
            fq_out = self.mgr.tmp_dir / self.mgr.append_filename(fq_in, "trimmed").name
            cmds = ["java", "-jar", "trimmomatic-0.40.jar", "SE", fq_in, fq_out, "ILLUMINA?", "..."]
            self.mgr.update_sample_name(i, fq_out)
            execute(cmds, f"Trimming reads for {fq_in}.")
        logger.info("Reads trimmed successfully!")
    
    def build_kallisto_index(self) -> Path:
        idx_file = f"{strip_all_extensions(self.mgr.reference_file)}.idx"
        cmds = ["kallisto", "index", "-i",
                self.mgr.kallisto_dir / idx_file, "-T", self.mgr.tmp_dir, "-t", self.mgr.p, self.mgr.reference_file]
        logger.debug(f"{idx_file}, {self.mgr.kallisto_dir}")
        logger.debug(cmds)

        execute(cmds, "Building kallisto index file.")
        logger.info("Kallisto index file built successfully!")
        return idx_file

    def kallisto_quantify(self, idx_file: Path):
        for item in os.listdir(self.mgr.kallisto_dir):
            if os.path.splitext(item) != ".idx":
                shutil.rmtree(self.mgr.kallisto_dir / item, ignore_errors=True)
        for sample, sample_name in self.mgr.sample_names.items():
            cmds = ["kallisto", "quant", "-i", self.mgr.kallisto_dir / idx_file, "-t", self.mgr.p, "--single", "-l", self.mgr.mean_frag_length,
                    "-s", self.mgr.std_frag_length, "-o", self.mgr.kallisto_dir / f"{sample_name}", sample]
            execute(cmds, f"Quantifying {sample_name} with kallisto.")
        logger.info("Samples quantified successfully!")

    def run_kallisto(self, idx_file: Path = None):
        if idx_file is None:
            idx_file = self.build_kallisto_index()
        self.kallisto_quantify(idx_file)

    def scale_abundances(self):
        cmds = ["Rscript", "deseq.R", self.mgr.annotation_file, self.mgr.alpha, self.mgr.l2FC_thresh, *self.mgr.samples]
        execute(cmds, "Tximporting...")