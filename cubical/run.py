import sys
import json
import time
import pathlib
import numpy as np
from math import comb

from percolation import cell_percolation

DEFAULT_DIM = 4
DEFAULT_TRIALS = 5000
DEFAULT_SCALES = {2: [10, 20, 40, 80], 4: [3, 4, 6, 8, 9]}


def run_scale(dim, L, homology, trials, seed=None):
    if seed is None:
        seed = 1000 * dim + L
    t0 = time.time()
    perc, giants, occ, ncells = cell_percolation(dim, L, trials, homology=homology,
                                                 seed=seed, verify=True)
    out = pathlib.Path(f"output/statistics/d{dim}_L{L}")
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "percentages", perc)
    np.save(out / "giants", giants)
    np.save(out / "occupation", occ)
    meta = {
        "dimension": dim, "homology": homology, "scale": L,
        "rank": comb(dim, homology), "iterations": trials, "cells": ncells,
        "model": "Bernoulli cell percolation",
        "complex": "Cubical", "seconds": round(time.time() - t0, 1),
    }
    json.dump(meta, open(out / "metadata.json", "w"), indent=2)
    print(f"L={L}  birth {perc[perc>0].mean():.4f}  giants {giants.mean():.2f}  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    dim, trials, homology, scales = DEFAULT_DIM, DEFAULT_TRIALS, None, []
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--dim":
            dim = int(args[i + 1]); i += 2
        elif a == "--homology":
            homology = int(args[i + 1]); i += 2
        elif a == "--trials":
            trials = int(args[i + 1]); i += 2
        else:
            scales.append(int(a)); i += 1
    if homology is None:
        homology = dim // 2
    if not scales:
        scales = DEFAULT_SCALES.get(dim, [3, 4, 6, 8, 9])
    for L in scales:
        run_scale(dim, L, homology, trials)
