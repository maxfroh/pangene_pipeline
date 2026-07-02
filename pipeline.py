#!/usr/bin/env python3
import argparse
import tomllib

import seaborn as sns

from src.python.pipeline_manager import PipelineManager


def main():
    sns.set_theme()
    open("out.txt", mode="w").close()
    open("err.txt", mode="w").close()

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", "-c", required=True, type=str, help="The config file to use."
    )
    parser.add_argument("--pangene", action="store_true")
    parser.add_argument("--benchmark", action="store_true")

    args = parser.parse_args()

    config_file = args.config

    with open(config_file, mode="rb") as cf:
        config_dict = tomllib.load(cf)

    # print(config_dict)

    pipeline = PipelineManager(config_dict)
    pipeline.setup()
    pipeline.run()


if __name__ == "__main__":
    main()
