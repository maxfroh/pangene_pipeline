import subprocess
import logging
import sys
import os

from typing import Any
from pathlib import Path


# logger for program (maybe temp)
logger = logging.Logger("Pipeline", level=logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s")
formatter.datefmt = "%Y-%m-%d %H:%M:%S"
handler.setFormatter(formatter)
logger.addHandler(handler)


def execute(cmds: list[str | Any], msg: str = None, out_file: str = "out.txt", err_file: str = "err.txt"):
    """
    Execute a line of code in the command line.

    :param cmds: The command to executed, passed as a list of strings. Non-strings will be converted before use with the builtin `str()` function.
    :type cmds: list[str | Any]
    :param msg: A message to log when the command is executed (optional).
    :type msg: str
    :param out_file: The file to redirect anything sent to stdout.
    :type out_file: str
    :param err_file: The file to redirect anything sent to stderr.
    :type err_file: str
    :raise: Any error that occurs during the command line execution of `cmds`.
    """

    if msg is not None:
        logger.info(msg)
    # must be strs
    cmds = [str(cmd) for cmd in cmds]
    logger.debug(f"Running: > {" ".join(cmds)}")
    with open(out_file, mode="a") as fout, open(err_file, mode="a") as errout:
        result = subprocess.run(
            cmds, stdout=fout, stderr=errout, check=True, text=True)
        try:
            result.check_returncode()
        except Exception as e:
            logger.error(f"Something went wrong! \n{e}")
            raise e


def gunzip(file: str):
    """
    Gunzips the provided file.

    :param file: A `.gz` file.
    :type file: str
    """

    cmds = ["gunzip", file]
    execute(cmds, f"Unzipping {file}.")


def strip_all_extensions(file: str | Path):
    """Remove all file extensions to just get the true basename of the path"""
    f = Path(file)
    # remove file extensions to just get name of sample
    return os.path.splitext(os.path.splitext(f.name)[0])[0]
