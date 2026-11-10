import itertools as it
import math
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure


def _fig_ax_if_neccessary(
    ax,
    subplot_kw: dict | None = None,
):
    if ax is None:
        return plt.subplots(constrained_layout=True, subplot_kw=subplot_kw)

    return ax.get_figure(), ax


def grid_layout(
    n_groups: int,
    *,
    n_cols: int = 2,
    col_width: float = 3,
    row_height: float = 4,
    sharex: bool = True,
    sharey: bool = False,
):
    n_cols = min(n_groups, n_cols)
    n_rows = math.ceil(n_groups / n_cols)
    fig, axs = plt.subplots(
        n_rows,
        n_cols,
        figsize=(n_cols * col_width, n_rows * row_height),
        sharex=sharex,
        sharey=sharey,
        layout="constrained",
        squeeze=False,
    )
    return fig, list(axs.flatten())


def plot_groups(
    groups: list[pd.DataFrame] | list[pd.Series],
    *,
    plot_kind: Literal[
        "line",
        "bar",
        "barh",
        "hist",
        "box",
        "kde",
        "density",
        "area",
        "pie",
        "scatter",
        "hexbin",
    ] = "line",
    n_cols: int = 2,
    col_width: float = 3,
    row_height: float = 4,
    sharex: bool = True,
    sharey: bool = False,
    grid: bool = True,
):
    fig, axs = grid_layout(
        len(groups),
        n_cols=n_cols,
        col_width=col_width,
        row_height=row_height,
        sharex=sharex,
        sharey=sharey,
    )

    for group, ax in zip(groups, axs):
        group.plot(ax=ax, kind=plot_kind, grid=grid)

    for i in range(len(groups), len(axs)):
        axs[i].set_visible(False)

    return fig, axs


###############################################################################
# Autogrid
###############################################################################


def _partition_by_order_of_magnitude(s: pd.Series) -> list[list[str]]:
    return [
        i.to_list()
        for i in np.floor(np.log10(s)).to_frame(name=0).groupby(0)[0].groups.values()  # type: ignore
    ]


def _split_large_groups(groups, max_size: int):
    return list(
        it.chain(
            *(
                (
                    [group]
                    if len(group) < max_size
                    else [  # type: ignore
                        list(i)
                        for i in np.array_split(group, math.ceil(len(group) / max_size))  # type: ignore
                    ]
                )
                for group in groups
            )
        )
    )  # type: ignore


def plot_autogrid(
    s: pd.Series | pd.DataFrame,
    *,
    plot_kind: Literal[
        "line",
        "bar",
        "barh",
        "hist",
        "box",
        "kde",
        "density",
        "area",
        "pie",
        "scatter",
        "hexbin",
    ] = "line",
    n_cols: int = 2,
    col_width: float = 4,
    row_height: float = 3,
    max_group_size: int = 6,
    grid: bool = True,
):
    if isinstance(s, pd.Series):
        group_names = _partition_by_order_of_magnitude(s)
    else:
        group_names = _partition_by_order_of_magnitude(s.max())

    group_names = _split_large_groups(group_names, max_size=max_group_size)

    groups: list[pd.Series] | list[pd.DataFrame]

    if isinstance(s, pd.Series):
        groups = [s.loc[group] for group in group_names]
    else:
        groups = [s.loc[:, group] for group in group_names]

    return plot_groups(
        groups,
        plot_kind=plot_kind,
        n_cols=n_cols,
        col_width=col_width,
        row_height=row_height,
        grid=grid,
    )


def savefig(fig: Figure, name: str) -> None:
    (folder := Path("out")).mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / name, bbox_inches="tight", dpi=300)
