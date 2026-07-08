#!/usr/bin/env python3
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from upsetplot import UpSet, util

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


# OVERRIDE!
class FixedUpSet(UpSet):
    def plot_matrix(self, ax):
        """Plot the matrix of intersection indicators onto ax"""
        ax = self._reorient(ax)
        data = self.intersections
        n_cats = data.index.nlevels

        inclusion = data.index.to_frame().values

        # Prepare styling
        styles = [
            [
                (
                    self.subset_styles[i]
                    if inclusion[i, j]
                    else {"facecolor": self._other_dots_color, "linewidth": 0}
                )
                for j in range(n_cats)
            ]
            for i in range(len(data))
        ]
        styles = sum(styles, [])  # flatten nested list
        style_columns = {
            "facecolor": "facecolors",
            "edgecolor": "edgecolors",
            "linewidth": "linewidths",
            "linestyle": "linestyles",
            "hatch": "hatch",
        }
        styles = (
            pd.DataFrame(styles)
            .reindex(columns=style_columns.keys())
            .astype(
                {
                    "facecolor": "O",
                    "edgecolor": "O",
                    "linewidth": float,
                    "linestyle": "O",
                    "hatch": "O",
                }
            )
        )
        styles["linewidth"] = styles["linewidth"].fillna(1)
        styles["facecolor"] = styles["facecolor"].fillna(self._facecolor)
        styles["edgecolor"] = styles["edgecolor"].fillna(styles["facecolor"])
        styles["linestyle"] = styles["linestyle"].fillna("solid")
        del styles["hatch"]  # not supported in matrix (currently)

        x = np.repeat(np.arange(len(data)), n_cats)
        y = np.tile(np.arange(n_cats), len(data))

        # Plot dots
        if self._element_size is not None:  # noqa
            s = (self._element_size * 0.35) ** 2
        else:
            # TODO: make s relative to colw
            s = 200
        ax.scatter(
            *self._swapaxes(x, y),
            s=s,
            zorder=10,
            **styles.rename(columns=style_columns),
        )

        # Plot lines
        if self._with_lines:
            idx = np.flatnonzero(inclusion)
            line_data = (
                pd.Series(y[idx], index=x[idx])
                .groupby(level=0)
                .aggregate(["min", "max"])
            )
            colors = pd.Series(
                [
                    style.get("edgecolor", style.get("facecolor", self._facecolor))
                    for style in self.subset_styles
                ],
                name="color",
            )
            line_data = line_data.join(colors)
            ax.vlines(
                line_data.index.values,
                line_data["min"],
                line_data["max"],
                lw=2,
                colors=line_data["color"],
                zorder=5,
            )

        # Ticks and axes
        tick_axis = ax.yaxis
        tick_axis.set_ticks(np.arange(n_cats))
        tick_axis.set_ticklabels(
            data.index.names, rotation=0 if self._horizontal else -90
        )
        ax.xaxis.set_visible(False)
        ax.tick_params(axis="both", which="both", length=0)
        if not self._horizontal:
            ax.yaxis.set_ticks_position("top")
        ax.set_frame_on(False)
        ax.set_xlim(-0.5, x[-1] + 0.5)
        ax.set_autoscale_on(False)
        ax.grid(False)

    def _label_sizes(self, ax, rects, where):
        if not self._show_counts and not self._show_percentages:
            return
        if self._show_counts is True:
            count_fmt = "{:.0f}"
        else:
            count_fmt = self._show_counts
            if "{" not in count_fmt:
                count_fmt = util.to_new_pos_format(count_fmt)

        pct_fmt = "{:.1%}" if self._show_percentages is True else self._show_percentages

        if count_fmt and pct_fmt:
            if where == "top":
                fmt = f"{count_fmt}\n({pct_fmt})"
            else:
                fmt = f"{count_fmt} ({pct_fmt})"

            def make_args(val):
                return val, val / self.total

        elif count_fmt:
            fmt = count_fmt

            def make_args(val):
                return (val,)

        else:
            fmt = pct_fmt

            def make_args(val):
                return (val / self.total,)

        if where == "right":
            margin = 0.01 * abs(np.diff(ax.get_xlim())[0])
            for rect in rects:
                width = rect.get_width() + rect.get_x()
                ax.text(
                    width + margin,
                    rect.get_y() + rect.get_height() * 0.5,
                    fmt.format(*make_args(width)),
                    ha="left",
                    va="center",
                )
        elif where == "left":
            margin = 0.01 * abs(np.diff(ax.get_xlim())[0])
            for rect in rects:
                width = rect.get_width() + rect.get_x()
                ax.text(
                    width + margin,
                    rect.get_y() + rect.get_height() * 0.5,
                    fmt.format(*make_args(width)),
                    ha="right",
                    va="center",
                )
        elif where == "top":
            margin = 0.01 * abs(np.diff(ax.get_ylim())[0])
            for rect in rects:
                height = rect.get_height() + rect.get_y()
                ax.text(
                    rect.get_x() + rect.get_width() * 0.5,
                    height + margin,
                    fmt.format(*make_args(height)),
                    ha="center",
                    va="bottom",
                )
        else:
            raise NotImplementedError("unhandled where: %r" % where)
