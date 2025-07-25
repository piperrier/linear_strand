from sage.matrix.special import zero_matrix, block_matrix
from sage.modules.free_module_element import vector
from copy import copy

#################################
### construction the expanded matrix
#################################


def construct_B_from_A(A_matrix):
    ring = A_matrix.base_ring()
    field = ring.base()
    nu = len(ring.gens())
    n_forme = A_matrix.ncols()
    n_base = A_matrix.nrows()
    B = zero_matrix(field, n_forme, n_base*nu)

    for i in range(0,n_forme):
        row = []
        for p in A_matrix[:,i]:
            vec = vector([0]*nu)
            for k,v in p[0].dict().items():
                vec+=vector(k)*v
            #row[:0] = vec
            row += list(vec)
        B[i,:] = vector(row)
    B.subdivide([],[nu*i for i in range(1,n_forme-2 )])
    return(B)


#################################
### Compute syzygetic matrix
#################################


def comp_block_ijk(matrix, nu, m, n, i, j, k):
    #i,j,k same notations as Albano & La Scala
    # 1 <= i <= m
    # 1 <= j <= n
    # 1 <= k <= nu
    field = matrix.base_ring()

    block_ijk = zero_matrix(field, nu-k, nu)

    for row_index in range(nu-k):
        block_ijk[row_index, row_index+k] = matrix[j-1,(i-1)*nu + k]
        block_ijk[row_index, k] = matrix[j-1, (i-1)*nu + k+row_index]

    return block_ijk


def comp_block_ij(matrix, nu, m, n, i, j):
    block_ij = comp_block_ijk(matrix, nu, m, n, i, j, 0)

    for k in range(1, nu):
        block_ijk = comp_block_ijk(matrix, nu, m, n, i, j, k)
        block_ij = block_ij.stack(block_ijk)

    return block_ij


def matrix_syzygetic(matrix, nu, m, n):
    k = matrix.base_ring()

    blocks = []
    for i in range(1, m+1):
        for j in range(1, n+1):
            blocks.append(comp_block_ij(matrix, nu, m, n, i, j))

    m_syz = block_matrix(blocks, nrows = m, ncols = n) #sparse=true
    return m_syz


#################################
### matrix reduction
#################################


def block_ideal_elim(matrix, r, k_ideal, t_details):
    """
    elim non zero elements in "ideal" blocks
    """
    subdivisions = matrix.subdivisions()[0]
    subdivisions.insert(0,0)

    # pour chaque bloque diag
    for i in range(len(subdivisions)-1):
        
        subsubdivision = t_details[i]
        subsubdivision.insert(0,0)
        row_end = subdivisions[i+1] #excluded

        #pour chaque pivot du bloque
        for j in range(len(subsubdivision)-1):
            row_begin = subdivisions[i] + subsubdivision[j]
            #for k in range(len(k_ideal[i])-1,j-1,-1):
            for k in range(j,len(k_ideal[i])):
            # FIXME: j might be false in the following formula
                col = subdivisions[i] + sum(subsubdivision[:j+1]) + k_ideal[i][k] - j
                for l in range(row_begin,row_end):
                    if l != col :
                        matrix[l,:] = matrix[l,:] - matrix[l,col]*matrix[col,:]
    return matrix


def block_forward_substitution(matrix):
    """
    elim every blocks below identity blocks by row operations
    """
    subdivisions = matrix.subdivisions()[0]
    subdivisions.insert(0,0)
    subdivisions.append(matrix.nrows())
    
    for i in range(len(subdivisions)-2):
        col_begin = subdivisions[i]
        col_end = subdivisions[i+1]
        for j in range(i+1, len(subdivisions)-1):
            row_begin = subdivisions[j]
            row_end = subdivisions[j+1]
            matrix[row_begin:row_end] = matrix[row_begin:row_end] - matrix[row_begin:row_end,col_begin:col_end]*matrix[col_begin:col_end,:]
    return matrix

#################################
### matrix subdivision
#################################

def t_subdivide(matrix,t):
    t_prime = [x for x in t if x!=0]
    subdivision = [sum(t_prime[:(i+1)]) for i in range(len(t_prime))]
    new_matrix = copy(matrix)
    new_matrix.subdivide(subdivision,subdivision)
    return new_matrix


def block_T_TT(matrix):
    subdivisions = matrix.subdivisions()[0]
    T = matrix[0:subdivisions[-1],subdivisions[-1]:]
    TT = matrix.subdivision(len(subdivisions), len(subdivisions))
    return(T,TT)