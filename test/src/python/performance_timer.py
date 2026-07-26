from time import perf_counter
from src.python.utils import build_logger
from pathlib import Path

# NB: turn off deg analysis (pipelinemanager.run?) 

class PerformanceTimer():
    def __init__(self):
        self.curr_run = None
        self.curr_ref = None
        self.curr_pangene = None
        self.logger = build_logger("Timer")
        self.time_dict = {"Pipeline": {}}    
        self.steps = []
    
    def set_curr_run(self, run: str):
        self.curr_run = run
        self.curr_ref = None
        if self.curr_run not in self.time_dict:
            self.time_dict[self.curr_run] = {}
            
    def set_curr_pangene(self, pangene: str):
        self.curr_pangene = pangene
        self.curr_run = None
        if self.curr_pangene not in self.time_dict:
            self.time_dict[self.curr_pangene] = {}
    
    def set_curr_ref(self, ref: str):
        self.curr_ref = ref
        if self.curr_ref not in self.time_dict:
            self.time_dict[self.curr_run][self.curr_ref] = {} 
    
    def set_to_pipeline(self):
        self.curr_ref = None
        self.curr_run = None
        self.curr_pangene = None
            
    def add_time(self, func_name: str, start: bool):
        time_point = "start" if start else "finish"
        time = perf_counter()
        if self.curr_run is None:
            if self.curr_pangene is None:
                self.time_dict["Pipeline"][func_name][time_point] = time
                self.steps.append(("Pipeline", func_name, time_point))
            else:
                self.time_dict[self.curr_pangene][func_name][time_point] = time
                self.steps.append((self.curr_pangene, func_name, time_point))
        else:
            if self.curr_ref is None:
                self.time_dict[self.curr_run][func_name][time_point] = time
                self.steps.append((self.curr_run, func_name, time_point))
            else:
                self.time_dict[self.curr_run][self.curr_ref][func_name][time_point] = time
                self.steps.append((self.curr_run, self.curr_ref, func_name, time_point))
    
    def checkpoint(self):
        self.logger.debug("Saving a time checkpoint!")
        self.write_log()
                
    def write_log(self, log_loc: str | Path = "time.log"):
        log_loc = Path(log_loc)
        with log_loc.open(mode="a") as f:
            for step in self.steps:
                data = None
                for loc in step:
                    if data is None:
                        data = self.time_dict[loc]
                    else:
                        data = data[loc]
                f.write(f"{' | '.join(step)}\t{data}\n")
        self.steps = []
                
                
PT = PerformanceTimer()