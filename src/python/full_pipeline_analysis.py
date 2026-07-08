#!/usr/bin/env python3
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
from matplotlib_set_diagrams import EulerDiagram

from .run_manager import RunManager
from .utils import build_logger


class FullPipelineAnalyzer:
    def __init__(self, runs_dict: dict[str, RunManager]):
        pass

    def analyze_runs(self):
        pass
