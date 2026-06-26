#!/usr/bin/env python3
import argparse
import tomllib

from src.defaults import _ALL_STEPS
from src.python.deg import DEG
from src.python.managers import ReferenceManager
from src.python.pangene_constructor import PangeneConstructor
from src.python.utils import *


def main(args: list):
    open("out.txt", mode="w").close()
    open("err.txt", mode="w").close()

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", "-c", required=True, type=str, help="The config file to use."
    )
    parser.add_argument("--pangene", action="store_true")
    parser.add_argument("--benchmark", action="store_true")

    args = parser.parse_args(args)

    config_file = args.config

    with open(config_file, mode="rb") as cf:
        config = tomllib.load(cf)

    print(config)

    managers = ReferenceManager.get_runs_and_references(config)
    for manager in managers:
        print(manager["run"], manager["reference"])

        # deg = DEG(manager)
        # deg.perform_deg()
        constructor = PangeneConstructor(manager)


if __name__ == "__main__":
    main()
