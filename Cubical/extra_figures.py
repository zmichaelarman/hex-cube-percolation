import pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import comb

from config import histograms
from cubical import build
from ateams.arithmetic import ComputePersistencePairs

plt.rcParams.update(**histograms.occupation.rcParams)
OUT = pathlib.Path("output/figures")
OUT.mkdir(parents=True, exist_ok=True)


def filtration_births(d, N, trials, homology, seed=0):
    C = build(d, N)
    rank = comb(d, homology)
    low, high = int(C.breaks[homology]), int(C.breaks[homology + 1])
    ncells = high - low
    cellCount = int(C.tranches[d][1])
    full = C.matrices.full
    times = set(range(cellCount))
    rng = np.random.default_rng(seed)
    perc = np.zeros((trials, rank))
    for t in range(trials):
        filt = np.arange(cellCount)
        filt[low:high] = low + rng.permutation(ncells)
        births, deaths = zip(*ComputePersistencePairs(full, filt, homology, C.breaks))
        ess = sorted(e for e in (times - (set(births) | set(deaths))) if low <= e < high)
        dens = sorted((e - low + 1) / ncells for e in ess)
        perc[t, :len(dens)] = dens[:rank]
    return perc


def transition_figure(dim, RANK, data):
    allb = np.concatenate([p[p > 0] for p, _ in data.values()])
    lo, hi = allb.min(), allb.max()
    ps = np.linspace(lo, hi, 400)
    Ls = sorted(data)
    colors = plt.cm.plasma(np.linspace(0, 0.8, len(Ls)))
    fig, ax = plt.subplots(figsize=histograms.occupation.figsize)
    for c, L in zip(colors, Ls):
        perc, tr = data[L]
        b = np.sort(perc[perc > 0])
        ax.plot(ps, np.searchsorted(b, ps, side="right") / tr, color=c, lw=1, label=rf"$L={L}$")
    ax.axvline(0.5, ls="--", lw=0.8, color="k", alpha=0.6)
    ax.axhline(RANK, ls=":", lw=0.8, color="k", alpha=0.4)
    ax.axhline(RANK / 2, ls=":", lw=0.6, color="k", alpha=0.25)
    ax.set_xlabel(r"Occupation probability $p$")
    ax.set_ylabel(r"Mean number of giant cycles")
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.2, RANK + 0.3)
    ax.legend(fontsize=7, loc="upper left", ncol=2)
    fig.savefig(OUT / f"transition.d{dim}.png", **histograms.occupation.savefig)
    plt.close(fig)


def scaling_figure(dim, data):
    Ls = np.array(sorted(data), dtype=float)
    width = np.array([data[L][0][data[L][0] > 0].std() for L in sorted(data)])
    b, a = np.polyfit(np.log(Ls), np.log(width), 1)
    nu = -1.0 / b
    sg = histograms.occupation.colors(3)[-1]
    fig, ax = plt.subplots(figsize=histograms.occupation.figsize)
    Lf = np.array([Ls.min(), Ls.max()])
    ax.plot(Lf, np.exp(a) * Lf ** b, ls="--", lw=1, color=sg, alpha=0.8,
            label=rf"$\sim L^{{{b:.2f}}}$ \ ($1/\nu={-b:.2f}$)")
    ax.scatter(Ls, width, marker="s", s=28, facecolor="w", edgecolors=sg,
               linewidths=1.1, zorder=10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(Ls)
    ax.set_xticklabels([f"${int(x)}$" for x in Ls])
    ax.set_xlabel(r"Lattice scale $L$")
    ax.set_ylabel(r"Birth-density width (std)")
    ax.legend(fontsize=8, loc="upper right")
    fig.savefig(OUT / f"scaling.d{dim}.png", **histograms.occupation.savefig)
    plt.close(fig)
    return b, nu


if __name__ == "__main__":
    DIM, H, SCALES, TRIALS = 4, 2, [3, 4, 6, 8, 9], 2000
    data = {}
    for L in SCALES:
        data[L] = (filtration_births(DIM, L, TRIALS, H, seed=L), TRIALS)
        print(f"L={L}  birth {data[L][0][data[L][0]>0].mean():.4f}", flush=True)
    RANK = comb(DIM, H)
    transition_figure(DIM, RANK, data)
    b, nu = scaling_figure(DIM, data)
    print(f"d={DIM}: width ~ L^{b:.2f}  (nu {nu:.2f})")
    print("figures written")
