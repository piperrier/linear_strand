# Linear strand computation

This repository contains an experimental implementation of the syzygy distinguisher presented in <https://eprint.iacr.org/2024/1193> (see <https://github.com/randriam/syzygies> for legacy magma code).

This is a toy implementation of the algorithm described in: [Albano & La Scala's paper](https://link.springer.com/article/10.1007/s002000000043)

## Description

We are interested in computing the linear strand of a graded modul over a polynomial ring until a given $r$.

Let $\mathcal{R}=\mathbb{K}[x_1 \ldots x_{\nu}]$ be a graded polynomial.

Let $\mathcal{F}$ be a graded free module over $\mathcal{R}$ and let $\{e_1 \ldots e_m\}$ be a free basis of $\mathcal{F}$, assume $deg(e_i)=d$.

Let $\mathcal{M}$ be a graded submodule of $\mathcal{F}$ and let $\{f_1 \ldots f_n\}$ be a minimal set of generators of $\mathcal{M}$, assume $deg(f_j)=d+1$.

Then $f_j = \sum_{i=1}^mL_{ij}e_i$ where $L_{ij}$ is a linear form. The $f_j$ 's are linearly independent over $\mathbb{K}$, we are interested in finding a basis $L_j'$ of all the linear-form-vectors such that : $\sum_{i=1}^nL_{ij}'f_i=0\,$ that is $A\times L_j' = 0$, with :

```math
A = \begin{bmatrix}
L_{11} & L_{12} & \ldots & L_{1n} \\
L_{21} & L_{22} & \ldots & L_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
L_{m1} & L_{m2} & \ldots & L_{mn} \\
\end{bmatrix}\in \mathcal{R}_1^{m \times n}

\hspace{1cm}

L_j' = \begin{bmatrix}
L_{1j}' \\
L_{2j}' \\
\vdots \\
L_{nj}' \\
\end{bmatrix} \in \mathcal{R}_1^{n \times 1}

\hspace{1cm}

B = \begin{bmatrix}
a_{11}^1 \ldots a_{11}^\nu & \ldots & a_{m1}^1 \ldots a_{m1}^\nu \\
a_{12}^1 \ldots a_{12}^\nu & \ldots & a_{m2}^1 \ldots a_{m2}^\nu \\
\vdots                     & \ddots & \vdots                     \\
a_{1n}^1 \ldots a_{1n}^\nu & \ldots & a_{mn}^1 \ldots a_{mn}^\nu \\
\end{bmatrix} \in \mathbb{K}^{n \times m\nu}
```

Such a thing can be done by constructing the macaulay matrix $M$ of the polynomial system and computing its kernel. In addition, we would like to exploit the sparsness and the structure of the macaulay matrix to construct a reduced matrix and thus improve kernel computation.

The procedure is as follows:

- First, construct at a minimal cost the matrices $T$ and $\bar{T}$ that correspond to the simplified system described by the macaulay matrix $M$:

```math
M_{red} = \begin{bmatrix}
Id_n & T \\
0 & \bar{T} \\
\end{bmatrix}
```

- Then, compute $Ker(\bar{T})$ and extend it to $Ker(S) = Vec \left< (-T\times x | x), x \in Ker(\bar{T}) \right>$ (mutliply by -T and concatenate) and then we can permute columns back to obtain $Ker(M)$

- Finally, iterate:
  - $m \leftarrow \\#(f_i)_i$
  - $(e_i)_i \leftarrow (f_i)_i$
  - $n \leftarrow \\#(L'_i)_i$
  - $(f_j)_j \leftarrow (\sum_{i=1}^nL_{ij}'f_i)_j$

## Some (minor) modifications of Albano & La Scala conventions

Notations are almost the same as the one from Albano & La Scala:

- We switched to the lexicographical ordering as it seems more natural and we didn't see any reason to keep the order Albano & La Scala introduced :  
   $(h, k) < (h' , k' )$ iff $(k > k')\lor\left((k = k') \land (h < h')\right)$ becomes $(h, k) < (h' , k' )$ iff $(h < h')\lor\left((h = h') \land (k < k')\right)$
- As a consequence, we adapt the writing of B, now we read B from left to right, and we demand B to be Left Row Echelon Form.

## Installation

You just need sage math.

## Usage

## Organisation & files

### `rs_10_5.ipynb`

We compute the linear strand of the (10, 5) Reed-Solomon code.  
The matrix A is the M23 macaulay matrix, the first matrix we iterate on to get further syzygies.  
In this example, we **construct the whole matrix** associated to the system and apply the procedure described in Albano & La Scala step by step.

### `code.ipynb`

A notebook to fill with examples, still incomplete and surprisingly slow.

### `syzygies.sage`

This is Nicolas's implementation.  
We **only construct the *solution columns*** of the system, transvections we can deduce from reading the matrix B are "done" (in fact we are not even doing them) directly when constructing the matrix of the system.

### `code.sage`, `combinatorics.sage`, `lin_betti.sage`, `matrix.sage`, `utils.sage`

This is Pierre's implementation.  
There might be too much, and many useless stuff.
Currently it contains a version that **construct the whole system** and the version that only construct *solution columns* is still in progress.

- `code.sage` : compute the square of the code and Macaulay23 matrix, **the function computing the Macaulay23 matrix is surprisingly slow** (see `code.ipynb`, really slow...).
- `combinatorics.sage` : monomial orders, rows and columns permutation
- `lin_betti.sage` : baby syzygies distinguisher, with Albano & La Scala method to compute linear strand, and some examples/tests.
- `matrix.sage` : matrix constructions, manipulations and operations
- `utils.sage` : usefull stuff, including preprocessing

## Potential optimizations

- parallelization/multiprocessing
- GPU
- dense/sparse algebra
- predictable execution flow
- predictable memory-access pattern
- improve row and column permutations
- characteristic of the field
- use `numpy` and `numba`

## References

- [La Scala's slides](https://www.math.rwth-aachen.de/~Viktor.Levandovskyy/filez/semcalg0910/lascala_resolution.pdf)
- [Albano & La Scala paper](https://link.springer.com/article/10.1007/s002000000043)
- [Eisenbud](https://www-users.cse.umn.edu/~reiner/REU/REU2019notes/2005_Book_TheGeometryOfSyzygies.pdf)
