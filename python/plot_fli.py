#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_fli_data(filename: Path, skiprows: int = 0):
    data = np.loadtxt(filename, comments="#", skiprows=skiprows)

    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError("A fájlnak legalább három oszlopot kell tartalmaznia: a, e, FLI.")

    a = data[:, 0]
    e = data[:, 1]
    fli = data[:, 2]

    finite = np.isfinite(fli)

    if not np.any(finite):
        raise ValueError("A fájl nem tartalmaz véges FLI értéket.")

    max_fli = np.max(fli[finite])

    # NaN és ±Inf helyettesítése
    fli[~finite] = max_fli

    return a, e, fli


def create_fli_grid(a, e, fli):
    a_values = np.unique(a)
    e_values = np.unique(e)

    grid = np.full((len(e_values), len(a_values)), np.nan)

    a_index = {v: i for i, v in enumerate(a_values)}
    e_index = {v: i for i, v in enumerate(e_values)}

    for ai, ei, fi in zip(a, e, fli):
        grid[e_index[ei], a_index[ai]] = fi

    return a_values, e_values, grid


def plot_map(a_values, e_values, grid, output):
    fig, ax = plt.subplots(figsize=(12, 4.5))

    im = ax.pcolormesh(
        a_values,
        e_values,
        grid,
        shading="auto",
        cmap="turbo",
        rasterized=True,
    )

    ax.set_xlabel(r"$a$ (AU)", fontsize=16)
    ax.set_ylabel(r"$e$", fontsize=16)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("FLI")

    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("-o", "--output", default="fli_map.png")
    parser.add_argument("--skiprows", type=int, default=0)

    args = parser.parse_args()

    a, e, fli = read_fli_data(Path(args.input), args.skiprows)
    a_values, e_values, grid = create_fli_grid(a, e, fli)
    plot_map(a_values, e_values, grid, args.output)


if __name__ == "__main__":
    main()
