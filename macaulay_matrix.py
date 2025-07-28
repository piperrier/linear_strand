from sage.matrix.special import zero_matrix, block_matrix
from sage.modules.free_module_element import vector

#################################
### construction the expanded matrix
#################################


def construct_B_from_A(A_matrix):
    """
    Compute the expanded matrix of a polynomial matrix. Both matrix are dense.
    """
    ring = A_matrix.base_ring()
    field = ring.base()
    nu = len(ring.gens())
    n = A_matrix.ncols()
    m = A_matrix.nrows()
    B = zero_matrix(field, n, m*nu)

    for i in range(0,n):
        row = []
        for p in A_matrix[:,i]:
            vec = vector([0]*nu)
            for k,v in p[0].dict().items():
                vec+=vector(k)*v
            row += list(vec)
        B[i,:] = vector(row)
    B.subdivide([],[nu*i for i in range(1,n-2 )])
    
    return(B)


#################################
### Compute syzygetic matrix
#################################


def comp_block_ijk(matrix, nu, m, n, i, j, k):
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
    """
    Compute the whole macaulay matrix in dense representation.
    """
    blocks = []

    for i in range(1, m+1):
        for j in range(1, n+1):
            blocks.append(comp_block_ij(matrix, nu, m, n, i, j))

    m_syz = block_matrix(blocks, nrows = m, ncols = n) #sparse=true
    
    return m_syz


#################################
### matrix reduction
#################################


def block_ideal_elim(matrix, k_ideal, spacing, subdiv):
    """
    elim non zero elements in "ideal" blocks
    """
    subdivisions = [0] + subdiv
    k_nz = [x for x in k_ideal if x!= []]

    # for each diagonal block
    for i in range(len(subdivisions)-1):
        subsubdivision = [0] + spacing[i]
        row_end = subdivisions[i+1] #excluded

        # for each pivots of the block (sub-block)
        for j in range(len(subsubdivision)-1):
            row_begin = subdivisions[i] + subsubdivision[j]

            # for each "inferior" pivots
            for k in range(j,len(k_nz[i])):
                col = subdivisions[i] + sum(subsubdivision[:j+1]) + k_nz[i][k] - j

                for l in range(row_begin,row_end):
                    if l != col :
                        matrix[l,:] -= matrix[l,col]*matrix[col,:]
    
    return matrix


def block_forward_substitution(matrix, subdiv):
    """
    elim every blocks below identity blocks by row operations
    """
    subdivisions =[0] + subdiv + [matrix.nrows()]
        
    for i in range(len(subdivisions)-2):
        col_begin = subdivisions[i]
        col_end = subdivisions[i+1]
        for j in range(i+1, len(subdivisions)-1):
            row_begin = subdivisions[j]
            row_end = subdivisions[j+1]
            matrix[row_begin:row_end] -= matrix[row_begin:row_end,col_begin:col_end]*matrix[col_begin:col_end,:]
    
    return matrix


#################################
### matrix subdivision
#################################


def t_subdivide(matrix,subdiv):
    """
    subdivide the matrix, essentially for pretty print
    """
    matrix.subdivide(subdiv,subdiv)

    return matrix


def block_T_TT(matrix,subdiv):
    T  = matrix[:subdiv[-1], subdiv[-1]:]
    TT = matrix[subdiv[-1]:, subdiv[-1]:]
    
    return(T,TT)