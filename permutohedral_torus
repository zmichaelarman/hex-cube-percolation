
from __future__ import annotations

import numpy as np
from scipy import sparse
from itertools import permutations, combinations, product
from math import comb



#  Lattice geometry in c-coordinates (integer combinations of v_1..v_d)

def _direction_deltas(d):
    """
    The d+1 'step' vectors used to walk along a Delaunay simplex, expressed in
    c-coordinates (the integer lattice Z^d in the basis v_1..v_d).

    Adding v_k (1 <= k <= d) moves c by +e_k.
    Adding v_{d+1} moves c by -(1,...,1)  since v_{d+1} = -(v_1+...+v_d).
    Returns a list of length d+1 of integer numpy vectors (directions 0..d).
    """
    deltas = []
    for k in range(d):                       # directions 0..d-1  ->  +e_k
        e = np.zeros(d, dtype=int)
        e[k] = 1
        deltas.append(e)
    deltas.append(-np.ones(d, dtype=int))    # direction d        ->  -(1..1)
    return deltas


def _top_delaunay_simplices(d, N):
    """
    Generate every top-dimensional Delaunay simplex of A*_d on the torus
    (Z/N)^d, as a set of frozensets of lattice points (each lattice point a
    tuple in (Z/N)^d).

    each top Delaunay simplex is:   { L, L+v_{pi(1)}, L+v_{pi(1)}+v_{pi(2)}, ..., L+v_{pi(1)}+...+v_{pi(d)} }
    for a base lattice point L and a permutation pi of {1,...,d+1}.  We reduce
    every vertex mod N and drop wrap-degenerate simplices (fewer than d+1
    distinct points), which only occur when N is too small.
    """
    deltas = _direction_deltas(d)
    simplices = set()
    n_degenerate = 0
    for base in product(range(N), repeat=d):
        L = np.array(base, dtype=int)
        for perm in permutations(range(d + 1)):
            verts = [tuple(L % N)]
            c = L.copy()
            for step in perm[:d]:            # first d of the d+1 directions
                c = c + deltas[step]
                verts.append(tuple(c % N))
            s = frozenset(verts)
            if len(s) == d + 1:
                simplices.add(s)
            else:
                n_degenerate += 1
    return simplices, n_degenerate


#  Build the full Delaunay complex and the dual permutohedral boundary maps

class PermutohedralTorus:
    """
    The permutohedral complex on T^d_N together with its boundary matrices.

    Attributes:
    d, N                : dimension and scale.
    delaunay[m]         : sorted list of Delaunay m-simplices (each a sorted
                          tuple of m+1 lattice points).  A Delaunay m-simplex is
                          the dual of a permutohedral (d-m)-cell.
    delaunay_index[m]   : { simplex : row/col index } for dimension m.
    perm_cells[j]       : sorted list of permutohedral j-cells.  Identical data
                          to delaunay[d-j] (a permutohedral j-cell *is* the set
                          of lattice points whose permutohedra meet in it).
    boundary[j]         : scipy CSR matrix d^V_j : C_j -> C_{j-1}
                          (rows = (j-1)-cells, cols = j-cells), entries in {-1,0,1}.
    """

    def __init__(self, d, N):
        self.d = d
        self.N = N

        #all top Delaunay simplices (dual to permutohedral vertices)
        top, n_deg = _top_delaunay_simplices(d, N)
        self.n_degenerate_top = n_deg

        # all Delaunay simplices of every dimension = all subsets, deduped
        delaunay = {m: set() for m in range(d + 1)}
        for s in top:
            pts = sorted(s)
            for size in range(1, d + 2):                # subset sizes 1..d+1
                for sub in combinations(pts, size):
                    delaunay[size - 1].add(tuple(sub))   # already sorted
        self.delaunay = {m: sorted(delaunay[m]) for m in range(d + 1)}
        self.delaunay_index = {
            m: {s: i for i, s in enumerate(self.delaunay[m])} for m in range(d + 1)
        }

        # Delaunay simplicial boundary maps  d^K_m : C_m -> C_{m-1}
        self._delaunay_boundary = {
            m: self._build_delaunay_boundary(m) for m in range(1, d + 1)
        }

        #permutohedral cells and boundary maps via duality / transpose
        #  permutohedral j-cell  <->  Delaunay (d-j)-simplex
        self.perm_cells = {j: self.delaunay[d - j] for j in range(d + 1)}
        #  d^V_j = ( d^K_{d-j+1} )^T  for j = 1..d
        self.boundary = {
            j: self._delaunay_boundary[d - j + 1].transpose().tocsr()
            for j in range(1, d + 1)
        }

    # Delaunay simplicial boundary (standard alternating-sign rule)
    def _build_delaunay_boundary(self, m):
        rows = len(self.delaunay[m - 1])
        cols = len(self.delaunay[m])
        ri, ci, data = [], [], []
        idx_lo = self.delaunay_index[m - 1]
        for simplex, j in self.delaunay_index[m].items():
            for i_del in range(m + 1):
                face = simplex[:i_del] + simplex[i_del + 1:]
                ri.append(idx_lo[face])
                ci.append(j)
                data.append((-1) ** i_del)
        return sparse.coo_matrix((data, (ri, ci)), shape=(rows, cols),
                                 dtype=np.int64).tocsr()

    #convenience
    def num_cells(self, j):
        return len(self.perm_cells[j])

    def cell_counts(self):
        return [self.num_cells(j) for j in range(self.d + 1)]

    def euler_characteristic(self):
        return sum((-1) ** j * self.num_cells(j) for j in range(self.d + 1))



#  Verification

def _rank_mod_p(M, p):
    """Rank of an integer matrix over the prime field F_p (dense)."""
    A = (np.asarray(M.todense()) % p).astype(np.int64)
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c] % p != 0:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)            # modular inverse (p prime)
        A[r] = (A[r] * inv) % p
        for i in range(rows):
            if i != r and A[i, c] % p != 0:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        r += 1
        if r == rows:
            break
    return r


def _rank_Q(M):
    if M.shape[0] == 0 or M.shape[1] == 0:
        return 0
    return int(np.linalg.matrix_rank(np.asarray(M.todense(), dtype=float)))


def betti_numbers(P, field="Q"):
    """Betti numbers of the permutohedral complex from its boundary matrices."""
    d = P.d
    rank = (lambda M: _rank_Q(M)) if field == "Q" else (lambda M: _rank_mod_p(M, field))
    rk = {j: rank(P.boundary[j]) for j in range(1, d + 1)}
    betti = []
    for j in range(d + 1):
        n = P.num_cells(j)
        r_j = rk.get(j, 0)                # rank d^V_j        (out of dim-j space)
        r_jp1 = rk.get(j + 1, 0)          # rank d^V_{j+1}
        betti.append(n - r_j - r_jp1)
    return betti


def boundary_squared_is_zero(P):
    for j in range(2, P.d + 1):
        prod = P.boundary[j - 1] @ P.boundary[j]
        if prod.nnz != 0 and not np.allclose(prod.data, 0):
            return False, j
    return True, None


def verify(P, verbose=True):
    d, N = P.d, P.N
    expected_betti = [comb(d, j) for j in range(d + 1)]
    bb_zero, bad = boundary_squared_is_zero(P)
    euler = P.euler_characteristic()
    bQ = betti_numbers(P, "Q")
    b2 = betti_numbers(P, 2)
    ok = (bb_zero and euler == 0 and bQ == expected_betti and b2 == expected_betti)
    if verbose:
        print(f"  cells by dim (perm 0..d): {P.cell_counts()}")
        print(f"  #top permutohedra (d-cells) = {P.num_cells(d)}   (expected N^d = {N**d})")
        print(f"  #permutohedral vertices (0-cells) = {P.num_cells(0)}   "
              f"(= #top Delaunay simplices)")
        print(f"  boundary^2 = 0 : {bb_zero}" + ("" if bb_zero else f"  (fails at j={bad})"))
        print(f"  Euler characteristic = {euler}  (torus expects 0)")
        print(f"  Betti /Q  = {bQ}")
        print(f"  Betti /F2 = {b2}")
        print(f"  expected  = {expected_betti}   (C(d,j))")
        print(f"  ===> {'PASS' if ok else 'FAIL'}")
    return ok


#  printing the boundary matrices

def print_complex(P, show_matrices=True, max_dim_for_full=None):
    d, N = P.d, P.N
    print("=" * 70)
    print(f"PERMUTOHEDRAL COMPLEX ON THE RHOMBIC TORUS   d={d}, N={N}")
    print(f"(A*_{d} / {N} A*_{d};  {N**d} permutohedra)")
    print("=" * 70)

    print("\nLattice basis  v_k = (1/(d+1)) 1 - e_k   in R^(d+1):")
    for k in range(d):
        v = np.full(d + 1, 1.0 / (d + 1)); v[k] -= 1.0
        print(f"  v_{k+1} = {np.round(v, 4)}")
    vlast = -np.array([(1.0/(d+1)) - (1.0 if i == k else 0.0)
                       for k in range(d) for i in [None]])  # placeholder
    vd1 = np.full(d + 1, 1.0 / (d + 1)); vd1[d] -= 1.0
    print(f"  v_{d+1} = {np.round(vd1, 4)}   = -(v_1+...+v_d)")

    print("\nPermutohedral cell counts:")
    names = {0: "vertices", 1: "edges", d - 1: "facets", d: "permutohedra"}
    for j in range(d + 1):
        tag = names.get(j, "")
        print(f"  dim {j:>2}: {P.num_cells(j):>6}   {tag}")

    if show_matrices:
        for j in range(1, d + 1):
            B = P.boundary[j]
            print(f"\nBoundary  d_{j} : C_{j} -> C_{j-1}   "
                  f"(shape {B.shape[0]} x {B.shape[1]}, nnz={B.nnz})")
            full = (max_dim_for_full is None) or (max(B.shape) <= max_dim_for_full)
            if full and max(B.shape) <= 40:
                print(_format_matrix(B))
            else:
                print("  (large; showing nonzero entries of first 6 columns)")
                Bc = B.tocsc()
                for col in range(min(6, B.shape[1])):
                    s, e = Bc.indptr[col], Bc.indptr[col + 1]
                    entries = {int(Bc.indices[t]): int(Bc.data[t]) for t in range(s, e)}
                    print(f"    col {col}: {entries}")


def _format_matrix(B):
    A = np.asarray(B.todense(), dtype=int)
    rows, cols = A.shape
    out = []
    header = "      " + "".join(f"{c:>3}" for c in range(cols))
    out.append(header)
    out.append("     +" + "-" * (3 * cols))
    for r in range(rows):
        line = f"{r:>4} |"
        for c in range(cols):
            v = A[r, c]
            line += "  ." if v == 0 else f"{v:>3}"
        out.append(line)
    return "\n".join(out)


#  Driver: find minimum working N

if __name__ == "__main__":
    print("Scanning for the minimum torus size N that yields the correct")
    print("topology (Betti = C(d,j)) in each dimension.\n")

    for d in (2, 3, 4):
        print(f"--- d = {d}  (expected Betti {[comb(d, j) for j in range(d+1)]}) ---")
        Ns = range(2, 7) if d == 2 else (range(2, 6) if d == 3 else range(2, 4))
        first_ok = None
        for N in Ns:
            try:
                P = PermutohedralTorus(d, N)
            except Exception as exc:               
                print(f"  N={N}: ERROR {exc!r}")
                continue
            bQ = betti_numbers(P, "Q")
            exp = [comb(d, j) for j in range(d + 1)]
            tag = "ok" if bQ == exp else "wrong"
            print(f"  N={N}: cells={P.cell_counts()}  Betti/Q={bQ}  [{tag}]")
            if bQ == exp and first_ok is None:
                first_ok = N
        print(f"  minimum working N for d={d}: {first_ok}\n")
