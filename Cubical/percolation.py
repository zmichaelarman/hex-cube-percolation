import numpy as np
from math import comb

from ateams.models import Bernoulli
from ateams import Chain
from ateams.arithmetic import ComputePersistencePairs

from cubical import build


def _verify_full(C, homology):
    d = C.dimension
    rank = comb(d, homology)
    low, high = int(C.breaks[homology]), int(C.breaks[homology + 1])
    cellCount = int(C.tranches[d][1])
    times = set(range(cellCount))
    filt = np.arange(cellCount)
    births, deaths = zip(*ComputePersistencePairs(C.matrices.full, filt, homology, C.breaks))
    ess = [e for e in (times - (set(births) | set(deaths))) if low <= e < high]
    assert len(ess) == rank, f"verify failed: H_{homology}={len(ess)}, expected {rank}"


def cell_percolation(d, N, trials, homology=None, seed=None, verify=True):
    if homology is None:
        homology = d // 2
    C = build(d, N)
    rank = comb(d, homology)
    low = int(C.breaks[homology])
    ncells = int(C.breaks[homology + 1]) - low
    if verify:
        _verify_full(C, homology)

    model = Bernoulli(C, dimension=homology)
    if seed is not None:
        model.RNG = np.random.default_rng(seed)

    percentages = np.zeros((trials, rank))
    giants = np.zeros(trials, dtype=int)
    occupation = np.zeros(trials)
    for t, (occ, gpos) in enumerate(Chain(model, steps=trials)):
        gpos = np.asarray(gpos, dtype=float)
        k = min(gpos.shape[0], rank)
        occupation[t] = float(occ.sum()) / occ.shape[0]
        giants[t] = int(gpos.shape[0])
        if k:
            percentages[t, :k] = np.sort((gpos - low) / ncells)[:k]
    return percentages, giants, occupation, ncells
