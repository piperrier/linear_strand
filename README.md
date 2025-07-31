# Linear strand computation

This repository contains an experimental implementation of the syzygy distinguisher presented in <https://eprint.iacr.org/2024/1193> (see <https://github.com/randriam/syzygies> for legacy magma code).

This is a toy implementation of the algorithm described in: [Albano & La Scala's paper](https://link.springer.com/article/10.1007/s002000000043)

## Description

We are interested in computing the linear strand of a graded module over a polynomial ring until a given $r$.

Let $\mathcal{R}=\mathbb{K}[x_1 \ldots x_{\nu}]$ be a graded polynomial ring.

Let $\mathcal{F}$ be a graded free module over $\mathcal{R}$ and let $(e_1, \dots, e_m)$ be a free basis of $\mathcal{F}$, assume $\mbox{deg}(e_i)=d$.

Let $\mathcal{M}$ be a graded submodule of $\mathcal{F}$ and let $(f_1, \dots ,f_n)$ be a minimal set of generators of $\mathcal{M}$, assume $\mbox{deg}(f_j)=d+1$.

Then $f_j = \sum_{i=1}^mL_{ij}e_i$ where $L_{ij}$ is a linear form. The $f_j$'s can be written column wise in a matrix $A$ with linear-forms coefficients:

```math
A = \begin{bmatrix}
L_{11} & L_{12} & \ldots & L_{1n} \\
L_{21} & L_{22} & \ldots & L_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
L_{m1} & L_{m2} & \ldots & L_{mn}
\end{bmatrix}\in \mathcal{R}_1^{m \times n}
```

Alternatively, the $f_j$'s can be written row wise in an expanded matrix $B$ with coefficient in $\mathbb{K}$:

```math
B = \begin{bmatrix}
a_{11}^1 \ldots a_{11}^\nu & \ldots & a_{m1}^1 \ldots a_{m1}^\nu \\
a_{12}^1 \ldots a_{12}^\nu & \ldots & a_{m2}^1 \ldots a_{m2}^\nu \\
\vdots                     & \ddots & \vdots                     \\
a_{1n}^1 \ldots a_{1n}^\nu & \ldots & a_{mn}^1 \ldots a_{mn}^\nu
\end{bmatrix} \in \mathbb{K}^{n \times m\nu}

\hspace{0.5cm}

L_{ij} = \sum_{k=1}^{\nu}a_{ij}^kx_k
```

The $f_j$ 's are linearly independent over $\mathbb{K}$, we are interested in finding a basis $L_j'$ of all the linear-form-vectors such that : $\sum_{i=1}^nL_{ij}'f_i=0$ that is $A\times L_j' = 0$, with :

```math
L_j' = \begin{bmatrix}
L_{1j}' \\
L_{2j}' \\
\vdots \\
L_{nj}'
\end{bmatrix} \in \mathcal{R}_1^{n \times 1}
```

Such a thing can be done by constructing the macaulay matrix $M$ of the polynomial system and computing its kernel. In addition, we would like to exploit the sparsness and the structure of the macaulay matrix to construct a reduced matrix and thus improve kernel computation.

The procedure is as follows:

- First, construct at a minimal cost the matrices $T$ and $\bar{T}$ that correspond to the simplified system described by the macaulay matrix $M$:

```math
M_{red} = \begin{bmatrix}
Id_n & T \\
0 & \bar{T}
\end{bmatrix}
```

- Then, compute $Ker(\bar{T})$ and extend it to $Ker(M_{red}) = Vec \left< (-T\times x | x), x \in Ker(\bar{T}) \right>$ (multiply by -T and concatenate) and then permute columns back to obtain $Ker(M)$

- Finally, iterate:
  - $m \leftarrow \\#(f_i)_i$
  - $(e_i)_i \leftarrow (f_i)_i$
  - $n \leftarrow \\#(L'_i)_i$
  - $(f_j)\_j \leftarrow (\sum_{i=1}^nL_{ij}'f_i)_j$

## Some (minor) modifications of Albano & La Scala conventions

Notations are almost the same as the one from Albano & La Scala:

- We switched to the lexicographical ordering as it seems more natural and we didn't see any reason to keep the order Albano & La Scala introduced :  
  $(h, k) < (h' , k' )$ iff $(k > k')\lor\left((k = k') \land (h < h')\right)$ becomes $(h, k) < (h' , k' )$ iff $(h < h')\lor\left((h = h') \land (k < k')\right)$
- As a consequence, the notation of matrix B is adjusted. Now, B is read from left to right and is required to be in Left Row Echelon Form.

## Installation

You just need sage math.

## Usage

In a `sage` interpreter, run :

```python
  from linear_strand import *

  G = matrix(GF(3),
  [[1, 0, 0, 0, 0, 0, 2, 0, 1, 2, 1],
  [0, 1, 0, 0, 0, 0, 1, 2, 2, 2, 1],
  [0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1],
  [0, 0, 0, 1, 0, 0, 1, 1, 0, 2, 2],
  [0, 0, 0, 0, 1, 0, 2, 1, 2, 2, 0],
  [0, 0, 0, 0, 0, 1, 0, 2, 1, 2, 2]])
  golay = CodeStrand(G)
  dist = golay.compute(r=4)
  print(dist)
```

## Organisation & files

### `rs_10_5.ipynb`

We compute the linear strand of the (10, 5) Reed-Solomon code.  
The matrix A is the M23 macaulay matrix, the first matrix we iterate on to get further syzygies.  
In this example, we **construct the whole matrix** associated to the system and apply the procedure described in Albano & La Scala step by step.

### `syzygies.sage`

This is Nicolas's implementation.  
We **only construct the *solution columns*** of the system, transvections we can deduce from reading the matrix B are "done" (in fact we are not even doing them) directly when constructing the matrix of the system. This version can be converted in a sparse one quiet easily.

### `code_ideal.py`, `expanded_matrix.py`, `linear_strand.py`, `macaulay_matrix.py`, `order.py`, `projects.sage`, `utils.py`

This is Pierre's implementation.  
Currently it contains a dense version that **construct the whole system** and it should be adapted to construct only the *solution columns* and exploit sparsness.

- `code_ideal.py`: compute the square of the code and Macaulay23 matrix.
- `expanded_matrix.py`: compute the expanded matrix and run `kosz_syz_mats` method.
- `linear_strand.py`: Toy syzygy distinguisher, with Albano & La Scala method.
- `macaulay_matrix.py`: matrix constructions, operations and reduction.
- `order.py`: monomial orders, rows and columns permutation.
- `projects.py`: examples of computations.
- `utils.py`: usefull stuff.

## Potential optimizations

- parallelization/multiprocessing
- GPU
- dense/sparse algebra
- predictable execution flow
- predictable memory-access pattern
- characteristic of the field
- use `numpy` and `numba`

## References

- [La Scala's slides](https://www.math.rwth-aachen.de/~Viktor.Levandovskyy/filez/semcalg0910/lascala_resolution.pdf)
- [Albano & La Scala paper](https://link.springer.com/article/10.1007/s002000000043)
- [Eisenbud](https://www-users.cse.umn.edu/~reiner/REU/REU2019notes/2005_Book_TheGeometryOfSyzygies.pdf)
