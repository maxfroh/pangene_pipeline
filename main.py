import argparse
import tomllib

from src.python.deg import DEG
from src.python.utils import *
from src.python.managers import RunManager

def main():
    open("out.txt", mode="w").close()
    open("err.txt", mode="w").close()
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", required=True,
                        type=str, help="The config file to use.")
    args = parser.parse_args()
    
    config_file = args.config
    
    with open(config_file, mode="rb") as cf:
        config = tomllib.load(cf)
        
    print(config)
    
    managers = RunManager.get_runs_and_references(config)
    for manager in managers:
        print(manager["run"], manager["reference"])
        
        deg = DEG(manager)
        deg.perform_deg()
    # file_mgr = FileManager(config)
    # params_mgr = ParamsManager(config)
    
    # for run in config["runs"]:
    #     for reference in config["runs"][run]["references"]:
    #         print(run, reference)        
    #         subset_file_mgr = file_mgr.get_subset_manager(run, reference)
    #         print(subset_file_mgr)
    


if __name__ == "__main__":
    main()
