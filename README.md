# hex-cube-percolation

Resulting figures I made in PDF file

I opted against using qhull and HighVoronoi.jl because I don't know Julia and qhull is computational expensive and no periodic boundary support.



Cubical: ATEAMS PHAT | fixed-p snapshot

Permutodedral: GUHDI lower-star | sweep

## Use

### 2D
```bash
python run.py --dim 2 10 20 40     # site percolation on the hexagonal lattice (H1)
python figures.py                  
```

### 4D
```bash
python run.py 3 4 6 8 9            # dim 4 is the default
python figures.py
```

| flag | meaning | default |
|---|---|---|
| `--dim D` | lattice dimension | 4 |
| `--homology H` | which giant cycles | middle, `D//2` |
| `--trials N` | samples per scale | 5000 |
| `--orientation O` | `rhombic` or `square` torus | `rhombic` |

### Note
"square" only works in dimension 2,3,4 because I hardcoded the orthogonal directions of the torus


## Citations

ATEAMS:( https://github.com/apizzimenti/ATEAMS )

Adams, A., Baek, J., & Davis, M. A. (2010). Fast high‐dimensional filtering using the permutohedral lattice. Computer Graphics Forum, 29(2), 753–762. https://doi.org/10.1111/j.1467-8659.2009.01645.x

Baek, J., & Adams, A. (2009). Some useful properties of the permutohedral lattice for Gaussian filtering [Technical Report]. Stanford University.

Bobrowski, O., & Skraba, P. (2020). Homological percolation: The formation of giant k-cycles [Preprint]. arXiv. https://arxiv.org/abs/2005.14011

Duncan, P., Kahle, M., & Schweinhart, B. (2020). Homological percolation on a torus: Plaquettes and permutohedra [Preprint]. arXiv. https://arxiv.org/abs/2011.11903

Kim, M. (2023). Algebraic characterization of the Voronoi cell structure of the A_n lattice [Preprint]. arXiv. https://doi.org/10.48550/arxiv.2304.10186

Maria, C., Boissonnat, J.-D., Glisse, M., & Yvinec, M. (2014). The GUDHI library: Simplicial complexes and persistent homology. In Mathematical Software – ICMS 2014 (LNCS 8592, pp. 167–174). Springer.

Panov, T., & Tril, V. (2025). Permutohedral complex and complements of diagonal subspace arrangements [Preprint]. arXiv. https://doi.org/10.48550/arxiv.2504.04183


