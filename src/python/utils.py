#!/usr/bin/env python3
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

default_logger = logging.Logger("Pipeline", level=logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
formatter.datefmt = "%Y-%m-%d %H:%M:%S"
handler.setFormatter(formatter)
default_logger.addHandler(handler)


def build_logger(name: str):
    lgr = logging.Logger(name, level=logging.DEBUG)
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    formatter.datefmt = "%Y-%m-%d %H:%M:%S"
    handler.setFormatter(formatter)
    lgr.addHandler(handler)
    return lgr


def execute(
    cmds: list[str | Any],
    msg: str = None,
    stdout: str | int | Any = "out.txt",
    stderr: str | int = "err.txt",
    stdin: str = None,
    log_out: bool = True,
    log_err: bool = True,
    path_out: bool = False,
    file_out: bool = False,
    logger: logging.Logger = default_logger,
):
    """
    Execute a line of code in the command line.

    :param cmds: The command to executed, passed as a list of strings. Non-strings will be converted before use with the builtin `str()` function.
    :type cmds: list[str | Any]
    :param msg: A message to log when the command is executed (optional).
    :type msg: str
    :param stdout: The file or location to redirect anything sent to stdout.
    :type stdout: str | Path | int
    :param stderr: The file or location to redirect anything sent to stderr.
    :type stderr: str | Path | int
    :param log_out: If stdout is going to a log file (default is True).
    :type log_out: bool
    :param log_err: If stderr is going to a log file (default is True).
    :type log_err: bool
    :param path_out: If stdout is going to a file that must be opened (just a Path is passed)
    :type path_out: bool
    :param file_out: If stdout is going to a different file (i.e., output of command must be saved to a file). Defaults to `False`, since stdout by default goes to a log.
    :type file_out: bool
    :raise: Any error that occurs during the command line execution of `cmds`.
    """
    if msg is not None:
        logger.info(msg)
    # must be strs
    cmds = [str(cmd) for cmd in cmds]
    logger.debug(f"Running: > {" ".join(cmds)}")
    if sum([log_out, file_out, path_out]) > 1:
        raise RuntimeError("Cannot save stdout output to more than one location!")
    if log_out and not (file_out or path_out):
        fout = open(stdout, mode="a")
    elif path_out and not (log_out or file_out):
        fout = open(stdout, mode="w")
    elif file_out and not (path_out or log_out):
        fout = stdout
    else:
        fout = subprocess.DEVNULL
    if log_err:
        ferr = open(stderr, mode="a")
    else:
        ferr = subprocess.PIPE
    if stdin is None:
        result = subprocess.run(cmds, stdout=fout, stderr=ferr, check=True, text=True)
    else:
        result = subprocess.run(cmds, stdin=stdin, stdout=fout, stderr=ferr, check=True)
    try:
        result.check_returncode()
    except Exception as e:
        logger.error(f"Something went wrong! \n{e}")
        raise e
    finally:
        if log_out or path_out:
            fout.close()
        if log_err:
            ferr.close()


def execute_quiet(
    cmds: list[str], msg: str = None, stdin=None, stdout=subprocess.DEVNULL
):
    # must be strs
    cmds = [str(cmd) for cmd in cmds]
    if stdin is None:
        result = subprocess.run(
            cmds, stdout=stdout, stderr=subprocess.PIPE, check=True, text=True
        )
    else:
        result = subprocess.run(
            cmds,
            stdin=stdin,
            stdout=stdout,
            stderr=subprocess.PIPE,
            check=True,
            text=True,
        )
    try:
        result.check_returncode()
    except Exception as e:
        print(f"Something went wrong! \n{e}")
        raise e


def gunzip(file: str | Path, out_file: Path = None, replace: bool = True):
    """
    Gunzips the provided file.

    :param file: A `.gz` file.
    :type file: str
    """
    cmds = ["gunzip"]
    if not replace:
        cmds.append("-k")
    if out_file is not None:
        cmds.append("-c")
        cmds.append(file)
        execute(
            cmds, f"Unzipping {file}.", stdout=out_file, path_out=True, log_out=False
        )
    else:
        cmds.append(file)
        execute(cmds, f"Unzipping {file}.")


def strip_filename(file: str | Path):
    """Remove all file extensions to just get the true basename of the path"""
    f = Path(file)
    # remove file extensions to just get name of sample
    return os.path.splitext(os.path.splitext(f.name)[0])[0]


def get_name_ext_and_is_gzip(file: str | Path) -> tuple[str, str, bool]:
    f = Path(file)
    name, ext = os.path.splitext(f.name)
    name, ext2 = os.path.splitext(name)
    if ext == ".gz":
        return name, ext2, True
    else:
        return name, ext, False
