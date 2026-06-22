import numpy as np
from fractions import Fraction
from itertools import permutations, product
from collections import deque


def _direction_deltas(d):
    deltas = [np.eye(d, dtype=int)[k] for k in range(d)]
    deltas.append(-np.ones(d, dtype=int))
    return deltas


_SQUARE_GEN = {
    2: [[1, -1], [1, 1]],
    3: [[-1, -1, 0], [-1, 0, -1], [0, -1, -1]],
    4: [[1, -1, 0, 0], [1, 1, -1, -1], [0, 0, 1, -1], [1, 1, 1, 1]],
}


def square_repetitions(d, N):
    B = np.array(_SQUARE_GEN[d], dtype=float)
    Q = np.diag(B @ (np.eye(d) - 1.0 / (d + 1)) @ B.T)
    return np.maximum(2, np.rint(N * np.sqrt(d / ((d + 1) * Q))).astype(int))


def _adjugate(M):
    d = len(M)
    A = [[Fraction(int(M[i][j])) for j in range(d)] +
         [Fraction(int(i == k)) for k in range(d)] for i in range(d)]
    det = Fraction(1)
    for c in range(d):
        piv = next(r for r in range(c, d) if A[r][c] != 0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det *= A[c][c]
        f = Fraction(1) / A[c][c]
        A[c] = [x * f for x in A[c]]
        for r in range(d):
            if r != c and A[r][c] != 0:
                g = A[r][c]
                A[r] = [A[r][j] - g * A[c][j] for j in range(2 * d)]
    inv = [[A[i][d + j] for j in range(d)] for i in range(d)]
    deti = int(det)
    adj = np.array([[int(deti * inv[i][j]) for j in range(d)] for i in range(d)], dtype=np.int64)
    return adj, deti


def _square_reducer(M):
    M = np.array(M, dtype=np.int64)
    adj, det = _adjugate(M)
    if det < 0:
        adj, det = -adj, -det
    Mt, adjT = M.T.copy(), adj.T.copy()

    def reduce(P):
        return P - ((P @ adjT) // det) @ Mt

    return reduce, det


def _delaunay_rhombic(d, N):
    deltas = _direction_deltas(d)
    raw = set()
    for base in product(range(N), repeat=d):
        L = np.array(base, dtype=int)
        for perm in permutations(range(d + 1)):
            verts, c = [tuple(L % N)], L.copy()
            for step in perm[:d]:
                c = c + deltas[step]
                verts.append(tuple(c % N))
            s = frozenset(verts)
            if len(s) == d + 1:
                raw.add(s)
    return raw


def _delaunay_square(d, N):
    deltas = np.array([np.array(x, dtype=np.int64) for x in _direction_deltas(d)])
    m = square_repetitions(d, N)
    M = np.array(_SQUARE_GEN[d], dtype=np.int64).T * m[None, :]
    reduce, det = _square_reducer(M)
    steps = np.vstack([deltas, -deltas])
    start = tuple(int(x) for x in reduce(np.zeros((1, d), dtype=np.int64))[0])
    seen, dq = {start}, deque([np.array(start, dtype=np.int64)])
    while dq:
        p = dq.popleft()
        for q in reduce(p[None, :] + steps):
            qt = tuple(int(x) for x in q)
            if qt not in seen:
                seen.add(qt)
                dq.append(q)
    cosets = np.array(sorted(seen), dtype=np.int64)
    raw = set()
    for perm in permutations(range(d + 1)):
        cur, cols = cosets.copy(), [cosets]
        for step in perm[:d]:
            cur = cur + deltas[step]
            cols.append(reduce(cur))
        for pts in np.stack(cols, axis=1):
            s = frozenset(tuple(int(x) for x in v) for v in pts)
            if len(s) == d + 1:
                raw.add(s)
    return raw


def delaunay_torus(d, N, orientation="rhombic"):
    raw = _delaunay_square(d, N) if orientation == "square" else _delaunay_rhombic(d, N)
    vertices = sorted({p for s in raw for p in s})
    vidx = {p: i for i, p in enumerate(vertices)}
    simplices = sorted(tuple(sorted(vidx[p] for p in s)) for s in raw)
    return vertices, simplices, vidx
