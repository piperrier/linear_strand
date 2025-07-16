# Linear strand computation
Notations are almost the same as the one from Albano & La Scala.

We aim at constructing at a minimal cost the matrices T and TT that correspond to the simplified system described by the macaulay matrix M :
 |||
 |----|----|
 | Id | T |
 | 0  | TT |

 Then we compute the kernel of TT and extend it (mutliply by -T, concatenate and permute columns) to recover the kernel of M.

## Some (minor) modifications of Albano & La Scala conventions
- We changed to the lexicographical ordering as it seems more natural and we didn't see any reason to keep the order Albano & La Scala introduced :  
   (h, k) < (h' , k' ) if and only if k > k' or (k = k' and h < h' ) becomes  
   (h, k) < (h' , k' ) if and only if h < h' or (h = h' and k < k' )
- As a consequence, we adapt the writing of B, now we read B from left to right, and we demand B to be Left Row Echelon Form.
## Organisation & files
### rs_10_5.ipynb
We compute the linear strand of the (10, 5) Reed-Solomon code.  
The matrix A is the M23 macaulay matrix, the firs matrix we iterate on to get further syzygies.  
In this example, we **construct the whole matrix** associated to the system and apply the procedure described in Albano & La Scala step by step.
### code.ipynb
A notebook to fill with examples, still incomplete and surprisingly slow.
### syzygies.sage
This is Nicolas's implementation.  
We **only construct the *solution columns*** of the system, transvections we can deduce from reading the matrix B are "done" (in fact we are not even doing them) directly when constructing the matrix of the system.
### code.sage, combinatorics.sage, lin_betti.sage, matrix.sage, utils.sage
This is Pierre's implementation.  
There might be too much, and many useless stuff.
Currently it contains a version that **construct the whole system** and the version that only construct *solution columns* is still in progress.
- code.sage : compute the square of the code and Macaulay23 matrix, **the function computing the Macaulay23 matrix is surprisingly slow** (see code.ipynb, really slow...).
- combinatorics.sage : monomial orders, rows and columns permutation
- lin_betti.sage : baby syzygies distinguisher, with Albano & La Scala method to compute linear strand, and some examples/tests.
- matrix.sage : matrix constructions, manipulations and operations
- utils.sage : usefull stuff, including preprocessing
## Potential optimizations (speedups, memory)
- [ ] parallel/multiprocessing
- [ ] GPU
- [ ] dense/sparse algebra
- [x] predictable execution flow
- [ ] predictable memory-access pattern
- [ ] improve row and column permutations

## References
- [La Scala slides](https://www.math.rwth-aachen.de/~Viktor.Levandovskyy/filez/semcalg0910/lascala_resolution.pdf)
- [Albano & La Scala paper](https://link.springer.com/article/10.1007/s002000000043)
- [Eisenbud](https://www-users.cse.umn.edu/~reiner/REU/REU2019notes/2005_Book_TheGeometryOfSyzygies.pdf)