from ateams.complexes import Cubical


def build(d, N):
    return Cubical().fromCorners([N] * d, periodic=True)
