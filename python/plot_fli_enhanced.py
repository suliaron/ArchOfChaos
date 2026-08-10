#!/usr/bin/env python3
"""
FLI-atlasz készítése háromoszlopos adatfájlból.

Bemeneti oszlopok:
    1. a   dimenziótlan fél nagytengely
    2. e   excentricitás
    3. FLI vagy log10(FLI)

Fő funkciók:
    - a átszámítása AU-ba: a_AU = a * a_Jupiter
    - NaN és ±Inf FLI értékek helyettesítése a legnagyobb véges FLI-vel
    - nyers FLI / log10(FLI) automatikus, heurisztikus felismerése
    - hiányzó rácspontok interpolációja
    - kézzel megadható vmin és vmax
    - PNG, PDF és EPS kimenet
    - Times New Roman jellegű publikációs megjelenés
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_JUPITER_A_AU = 5.2044


def parse_arguments():
    """Parancssori argumentumok feldolgozása."""
    parser = argparse.ArgumentParser(
        description="FLI-atlasz készítése az (a,e) síkon."
    )

    parser.add_argument(
        "input_file",
        type=Path,
        help="A háromoszlopos bemeneti fájl: a, e, FLI.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("fli_map.png"),
        help="Elsődleges kimeneti fájl. Alapértelmezés: fli_map.png",
    )
    parser.add_argument(
        "--skiprows",
        type=int,
        default=0,
        help="A kihagyandó fejlécsorok száma.",
    )
    parser.add_argument(
        "--jupiter-a",
        type=float,
        default=DEFAULT_JUPITER_A_AU,
        help=(
            "A Jupiter fél nagytengelye AU-ban. "
            f"Alapértelmezés: {DEFAULT_JUPITER_A_AU}"
        ),
    )
    parser.add_argument(
        "--fli-scale",
        choices=("auto", "raw", "log"),
        default="auto",
        help=(
            "A 3. oszlop értelmezése: auto, raw vagy log. "
            "A log azt jelenti, hogy az oszlop már log10(FLI)."
        ),
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="A színskála alsó határa.",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="A színskála felső határa.",
    )
    parser.add_argument(
        "--no-interpolation",
        action="store_true",
        help="A hiányzó rácspontok interpolációjának kikapcsolása.",
    )
    parser.add_argument(
        "--save-all",
        action="store_true",
        help="Mentés PNG, PDF és EPS formátumban is.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Az ábra ne jelenjen meg interaktív ablakban.",
    )

    return parser.parse_args()


def read_data(filename: Path, skiprows: int, jupiter_a: float):
    """
    Beolvassa az adatokat, AU-ba konvertálja az a értékeket, és kezeli
    a NaN/Inf FLI értékeket.
    """
    if not filename.is_file():
        raise FileNotFoundError(f"A bemeneti fájl nem található: {filename}")

    data = np.loadtxt(filename, comments="#", skiprows=skiprows)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 3:
        raise ValueError(
            "A fájlnak legalább három oszlopot kell tartalmaznia: a, e, FLI."
        )

    a_dimensionless = data[:, 0]
    e = data[:, 1]
    fli = data[:, 2].astype(float, copy=True)

    if not np.all(np.isfinite(a_dimensionless)):
        raise ValueError("Az a oszlop NaN vagy Inf értéket tartalmaz.")

    if not np.all(np.isfinite(e)):
        raise ValueError("Az e oszlop NaN vagy Inf értéket tartalmaz.")

    finite_mask = np.isfinite(fli)

    if not np.any(finite_mask):
        raise ValueError("Nincs egyetlen véges FLI érték sem a fájlban.")

    largest_finite_fli = np.max(fli[finite_mask])
    replaced_count = np.count_nonzero(~finite_mask)
    fli[~finite_mask] = largest_finite_fli

    a_au = a_dimensionless * jupiter_a

    return a_au, e, fli, replaced_count, largest_finite_fli


def detect_fli_scale(fli: np.ndarray):
    """
    Heurisztikusan eldönti, hogy a bemenet nyers FLI vagy log10(FLI).

    Ez nem lehet minden adathalmaznál teljesen biztos, ezért az eredmény
    felülbírálható a --fli-scale raw vagy --fli-scale log kapcsolóval.
    """
    finite = fli[np.isfinite(fli)]

    minimum = np.min(finite)
    maximum = np.max(finite)
    median = np.median(finite)

    if minimum >= -100.0 and maximum <= 100.0 and median <= 50.0:
        return "log"

    return "raw"


def prepare_fli_for_plot(fli: np.ndarray, mode: str):
    """Előkészíti az FLI értékeket a log10(FLI) színskálához."""
    detected_mode = detect_fli_scale(fli) if mode == "auto" else mode

    if detected_mode == "log":
        return fli.copy(), detected_mode

    positive_mask = fli > 0.0

    if not np.any(positive_mask):
        raise ValueError(
            "A nyers FLI értékek között nincs pozitív érték, ezért "
            "a log10 transzformáció nem végezhető el."
        )

    smallest_positive = np.min(fli[positive_mask])
    adjusted = fli.copy()
    adjusted[~positive_mask] = smallest_positive

    return np.log10(adjusted), detected_mode


def create_grid(a, e, values):
    """A pontlistából kétdimenziós rácsot készít."""
    a_values = np.unique(a)
    e_values = np.unique(e)

    grid = np.full((len(e_values), len(a_values)), np.nan, dtype=float)

    a_index = {value: index for index, value in enumerate(a_values)}
    e_index = {value: index for index, value in enumerate(e_values)}

    for ai, ei, value in zip(a, e, values):
        ia = a_index[ai]
        ie = e_index[ei]

        if np.isfinite(grid[ie, ia]):
            raise ValueError(
                f"Az ({ai:.16g}, {ei:.16g}) ponthoz több érték tartozik."
            )

        grid[ie, ia] = value

    return a_values, e_values, grid


def interpolate_missing(grid: np.ndarray):
    """
    A hiányzó rácspontokat sor- és oszlopirányú lineáris interpolációval tölti ki.
    """
    result = grid.copy()
    rows, cols = result.shape

    for i in range(rows):
        row = result[i]
        valid = np.isfinite(row)

        if np.count_nonzero(valid) >= 2:
            x = np.arange(cols)
            row[~valid] = np.interp(x[~valid], x[valid], row[valid])

    for j in range(cols):
        column = result[:, j]
        valid = np.isfinite(column)

        if np.count_nonzero(valid) >= 2:
            y = np.arange(rows)
            column[~valid] = np.interp(y[~valid], y[valid], column[valid])

    # Maradék, izolált pontok kitöltése közvetlen szomszédok átlagával.
    for _ in range(rows + cols):
        missing = np.argwhere(~np.isfinite(result))

        if missing.size == 0:
            break

        changed = False

        for i, j in missing:
            neighbours = []

            for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                ni = i + di
                nj = j + dj

                if 0 <= ni < rows and 0 <= nj < cols:
                    value = result[ni, nj]

                    if np.isfinite(value):
                        neighbours.append(value)

            if neighbours:
                result[i, j] = np.mean(neighbours)
                changed = True

        if not changed:
            break

    if np.any(~np.isfinite(result)):
        finite_values = result[np.isfinite(result)]

        if finite_values.size == 0:
            raise ValueError("Az egész FLI-rács hiányzik.")

        result[~np.isfinite(result)] = np.max(finite_values)

    return result


def configure_style():
    """Publikációs jellegű grafikai beállítások."""
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 14,
            "axes.labelsize": 21,
            "xtick.labelsize": 15,
            "ytick.labelsize": 15,
        }
    )


def save_figure(fig, output: Path, save_all: bool):
    """Az ábra mentése a kért formátumban."""
    output.parent.mkdir(parents=True, exist_ok=True)

    if save_all:
        stem = output.with_suffix("")

        for extension in (".png", ".pdf", ".eps"):
            target = stem.with_suffix(extension)
            kwargs = {"bbox_inches": "tight"}

            if extension == ".png":
                kwargs["dpi"] = 300

            fig.savefig(target, **kwargs)
            print(f"Elmentve: {target}")
    else:
        kwargs = {"bbox_inches": "tight"}

        if output.suffix.lower() == ".png":
            kwargs["dpi"] = 300

        fig.savefig(output, **kwargs)
        print(f"Elmentve: {output}")


def plot_map(
    a_values,
    e_values,
    grid,
    output,
    vmin,
    vmax,
    jupiter_a,
    save_all,
    show_plot,
):
    """Az FLI-térkép elkészítése."""
    configure_style()

#    fig, ax = plt.subplots(figsize=(12.0, 5.0))
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_box_aspect(1)

    image = ax.pcolormesh(
        a_values,
        e_values,
        grid,
        shading="auto",
        #cmap="turbo",
        cmap="YlGnBu_r",
        vmin=vmin,
#        vmax=vmax,
        vmax=30,
        rasterized=True,
    )

    # Jupiter excentricitása
    e_j = 0.0489
    
    q_j = jupiter_a * (1.0 - e_j)
    Q_j = jupiter_a * (1.0 + e_j)
    
    e_curve = np.linspace(e_values.min(), e_values.max(), 500)
    
    # q = qJ
    a_q = q_j / (1.0 - e_curve)
    mask = (a_q >= a_values.min()) & (a_q <= a_values.max())
    ax.plot(a_q[mask], e_curve[mask], color="black", lw=1.5, linestyle="-", label=r"$q_J$")
    #ax.plot(a_q[mask], e_curve[mask], color="limegreen", lw=1)
    
    # Q = QJ
    a_Q = Q_j / (1.0 + e_curve)
    mask = (a_Q >= a_values.min()) & (a_Q <= a_values.max())
    ax.plot(a_Q[mask], e_curve[mask], color="black", lw=1.5, linestyle="-", label=r"$Q_J$",)
    #ax.plot(a_Q[mask], e_curve[mask], color="limegreen", lw=1)
    
    # Tisserand T=3
    a_t = np.linspace(a_values.min(), a_values.max(), 1200)
    
    term = ((3.0 - jupiter_a / a_t) / 2.0) ** 2 * (jupiter_a / a_t)
    e2 = 1.0 - term
    mask = e2 >= 0.0
    ax.plot(a_t[mask], np.sqrt(e2[mask]), color="black", lw=1.5, linestyle="--", label=r"$T_J=3$")

    ax.set_xlabel(r"$a\;(\mathrm{AU})$")
    ax.set_ylabel(r"$e$", rotation=0, labelpad=18)
    #ax.set_xlim(a_values.min(), a_values.max())
    #ax.set_ylim(e_values.min(), e_values.max())
    ax.set_xlim(np.floor(a_values.min()),np.ceil(a_values.max()))
    ax.set_ylim(0.0, 1.0)
    ax.tick_params(direction="out", length=6, width=1)

    ax.legend(loc="lower right", fontsize=11, frameon=True, facecolor="white", edgecolor="black")
 
    colorbar = fig.colorbar(image, ax=ax, pad=0.02, aspect=30, extend="max")
    colorbar.ax.tick_params(labelsize=13)
    colorbar.set_label(r"$\ln(\mathrm{FLI})$", fontsize=18)

    vmax_real = np.nanmax(grid)
    colorbar.ax.text(0.5, 1.08, rf"$\max={vmax_real:.1f}$", ha="center", va="bottom", fontsize=12, transform=colorbar.ax.transAxes)

    fig.tight_layout()
    save_figure(fig, output, save_all)

    if show_plot:
        plt.show()
    else:
        plt.close(fig)


def main():
    """A program belépési pontja."""
    args = parse_arguments()

    if args.jupiter_a <= 0.0:
        raise ValueError("A Jupiter fél nagytengelyének pozitívnak kell lennie.")

    a, e, fli, replaced_count, replacement_value = read_data(
        args.input_file,
        args.skiprows,
        args.jupiter_a,
    )

#    plotted_fli, detected_scale = prepare_fli_for_plot(
#        fli,
#        args.fli_scale,
#    )

    a_values, e_values, grid = create_grid(a, e, fli)
    missing_count = np.count_nonzero(~np.isfinite(grid))

    if missing_count > 0 and not args.no_interpolation:
        grid = interpolate_missing(grid)

    print(f"Jupiter fél nagytengelye: {args.jupiter_a:.8f} AU")
    print(f"a-tartomány: {a_values.min():.8f}–{a_values.max():.8f} AU")
    print(f"e-tartomány: {e_values.min():.8f}–{e_values.max():.8f}")
#    print(f"FLI-skála: {detected_scale}")
    print(f"NaN/Inf helyettesítések száma: {replaced_count}")

    if replaced_count > 0:
        print(f"Helyettesítő FLI érték: {replacement_value:.8e}")

    print(f"Hiányzó rácspontok száma: {missing_count}")
    print(f"Ábrázolt minimum: {np.nanmin(grid):.8e}")
    print(f"Ábrázolt maximum: {np.nanmax(grid):.8e}")

    plot_map(
        a_values,
        e_values,
        grid,
        args.output,
        args.vmin,
        args.vmax,
        args.jupiter_a,
        args.save_all,
        not args.no_show,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"Hiba: {error}", file=sys.stderr)
        sys.exit(1)
