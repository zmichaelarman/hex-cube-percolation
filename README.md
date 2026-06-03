# permutohedral-complex

tests so far:
--- d = 2  (expected Betti [1, 2, 1]) ---
  N=2: cells=[4, 6, 4]  Betti/Q=[1, 0, 1]  [wrong]
  N=3: cells=[18, 27, 9]  Betti/Q=[1, 2, 1]  [ok]
  N=4: cells=[32, 48, 16]  Betti/Q=[1, 2, 1]  [ok]
  N=5: cells=[50, 75, 25]  Betti/Q=[1, 2, 1]  [ok]
  N=6: cells=[72, 108, 36]  Betti/Q=[1, 2, 1]  [ok]
  minimum working N for d=2: 3

--- d = 3  (expected Betti [1, 3, 3, 1]) ---
  N=2: cells=[24, 48, 28, 8]  Betti/Q=[0, 3, 0, 1]  [wrong]
  N=3: cells=[162, 324, 189, 27]  Betti/Q=[1, 3, 3, 1]  [ok]
  N=4: cells=[384, 768, 448, 64]  Betti/Q=[1, 3, 3, 1]  [ok]
  N=5: cells=[750, 1500, 875, 125]  Betti/Q=[1, 3, 3, 1]  [ok]
  minimum working N for d=3: 3

--- d = 4  (expected Betti [1, 4, 6, 4, 1]) ---
  N=2: cells=[192, 480, 400, 120, 16]  Betti/Q=[1, 0, 6, 0, 1]  [wrong]
  N=3: cells=[1944, 4860, 4050, 1215, 81]  Betti/Q=[1, 4, 6, 4, 1]  [ok]
  minimum working N for d=4: 3




Not ready to upload, to be officially implemented to ATEAMS:( https://github.com/apizzimenti/ATEAMS ) in near future.

A Python implementation for generating boundary matrices for permutohedral/hexagonal complexes in arbitrary dimensions using simplicial decomposition.

Work in Progress.

Adams, A., Baek, J., & Davis, M. A. (2010). Fast High‐Dimensional Filtering Using the Permutohedral Lattice. Computer Graphics Forum, 29, 753–762. https://doi.org/10.1111/j.1467-8659.2009.01645.x

Panov, T., & Tril, V. (2025). Permutohedral complex and complements of diagonal subspace arrangements. arXiv. https://doi.org/10.48550/arxiv.2504.04183

Kim, M. (2023). Algebraic Characterization of the Voronoi Cell Structure of the A_n Lattice. arXiv. https://doi.org/10.48550/arxiv.2304.10186
