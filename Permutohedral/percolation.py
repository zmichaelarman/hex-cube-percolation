import numpy as np
import gudhi as gd
from math import comb
from itertools import combinations

from permutohedral import delaunay_torus


def build_tree(d, N, homology, verify=True, orientation="rhombic"):
    vertices, simplices, _ = delaunay_torus(d, N, orientation)
    Nv = len(vertices)
    st = gd.SimplexTree()
    W = homology + 2
    for s in simplices:
        for face in combinations(s, W):
            st.insert(list(face))
    if verify:
        st.compute_persistence(persistence_dim_max=True)
        betti = st.betti_numbers()
        expected = comb(d, homology)
        assert len(betti) > homology and betti[homology] == expected, \
            f"verify failed: betti={betti}, expected H_{homology}={expected}"
    return st, Nv


def site_percolation(d, N, trials, homology=None, seed=None, verify=True, orientation="rhombic"):
    if homology is None:
        homology = d // 2
    rank = comb(d, homology)
    st, Nv = build_tree(d, N, homology, verify=verify, orientation=orientation)

    all_s = [list(s) for s, _ in st.get_simplices()]
    W = homology + 2
    S_pad = np.array([s + [s[0]] * (W - len(s)) for s in all_s], dtype=np.int64)

    rng = np.random.default_rng(seed)
    percentages = np.zeros((trials, rank))
    giants = np.zeros(trials, dtype=int)
    for t in range(trials):
        vrank = rng.permutation(Nv)
        vals = vrank[S_pad].max(axis=1).astype(float)
        for s, fv in zip(all_s, vals):
            st.assign_filtration(s, fv)
        st.compute_persistence(persistence_dim_max=True)
        ess = sorted((b + 1) / Nv for (b, dth) in
                     st.persistence_intervals_in_dimension(homology) if dth == float("inf"))
        k = min(len(ess), rank)
        percentages[t, :k] = ess[:k]
        giants[t] = int(np.sum(np.asarray(ess) <= 0.5))
    occupation = np.full(trials, 0.5)
    return percentages, giants, occupation, Nv
