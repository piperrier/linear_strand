# ============================
### construction
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

# ============================
### Compute syzygetic matrix
def comp_block_ijk(matrix, nu, m,n,i, j, k):
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

def comp_block_ij(matrix,nu,m,n,i,j):
    block_ij = comp_block_ijk(matrix,nu, m,n,i,j,0)

    for k in range(1,nu):
        block_ijk = comp_block_ijk(matrix,nu,m,n,i,j,k)
        block_ij = block_ij.stack(block_ijk)

    return block_ij

def matrix_syzygetic(matrix,nu,m,n):
    k = matrix.base_ring()

    blocks = []
    for i in range(1,m+1):
        for j in range(1,n+1):
            blocks.append(comp_block_ij(matrix,nu,m,n,i,j))

    m_syz = block_matrix(blocks, nrows = m, ncols=n) #sparse=true
    return m_syz

# ===============================
### x_h * l_ij kth block / hth column
def comp_col_ijk_h(matrix,nu,m,n,i,j,k,h):
    field = matrix.base_ring()
    col_ijk_h = zero_matrix(field, nu-k, 1)
    
    if h==k:
        for row_index in range(nu-k):
            col_ijk_h[row_index, 0] = matrix[j, i*nu + k + row_index]
        
    elif h>k:
        col_ijk_h[h-k, 0] = matrix[j,i*nu + k]

    return col_ijk_h

def comp_col_ij_h(matrix,nu,m,n,i,j,h):
    col_ij_h = comp_col_ijk_h(matrix,nu, m,n,i,j,0,h)

    for k in range(1,nu):
        col_ijk_h = comp_col_ijk_h(matrix,nu,m,n,i,j,k,h)
        col_ij_h = col_ij_h.stack(col_ijk_h)

    return col_ij_h

def comp_col_j_h(matrix,nu,m,n,j,h):
    col_j_h = comp_col_ij_h(matrix,nu,m,n,0,j,h)

    for i in range(1,m):
        col_ij_h = comp_col_ij_h(matrix,nu,m,n,i,j,h)
        col_j_h = col_j_h.stack(col_ij_h)
    
    return col_j_h

# ===================================
### comp T, TT
def comp_T_TT(matrix,nu,m,n,r,k_ideal,t):
    cols_sol = col_sol(r,k_ideal,nu)
    print(cols_sol)
    if len(cols_sol) == 0 :
        return 0
    else :
        j,h = cols_sol[0]
        M = comp_col_j_h(matrix,nu,m,n,j,h)
        for (j,h) in cols_sol[1:]:
            M = M.augment(comp_col_j_h(matrix,nu,m,n,j,h))
    
    new_row_order = koszul_order_block(nu,m,n,k_ideal)
    M = block_kosz_row_perm(M,nu,new_row_order)
    #print(M)
    #print("")
    T = M[:sum(t),:]
    TT = M[sum(t):,:]

    #row_perm
    return (T,TT)
    
# ============================
### elim non zero elements in "ideal" blocks
def block_ideal_elim(matrix, r, k_ideal, t_details):
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
    



# ============================
### elim every blocks below identity blocks by row operations
# this could be parellelized
def block_forward_substitution(matrix):
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

# ==============================
###
def t_subdivide(matrix,t):
    t_prime = [x for x in t if x!=0]
    subdivision = [sum(t_prime[:(i+1)]) for i in range(len(t_prime))]
    new_matrix = copy(matrix)
    new_matrix.subdivide(subdivision,subdivision)
    return new_matrix

# ==============================
### Get T and TT
def block_T_TT(matrix):
    subdivisions = matrix.subdivisions()[0]
    T = matrix[0:subdivisions[-1],subdivisions[-1]:]
    TT = matrix.subdivision(len(subdivisions), len(subdivisions))
    return(T,TT)

# =============================
def comp_Ki(matrix,nu,m,n,r,k_ideal,i):
    field = matrix.base_ring()
    r_i = r[i]
    k_i = k_ideal[i]
    size = r_i*(r_i-1)/2
    Ki = zero_matrix(field, r_i*nu - size, size)
    
    # forme lineaire
    for j in range(1,r_i):
        # col sol
        for l in range(j):
            shift = sum([x for x in range(j)])
            L2 = vector(matrix[sum(r[:i]) + j, i*nu + k_i[j]:(i+1)*nu])
            L1 = vector(-matrix[sum(r[:i])+ l, i*nu + k_i[j]:(i+1)*nu])            
            L1 = vector([coord for index, coord in enumerate(L1) if index not in k_ideal[:j]])

            #FIXME: -l might be false maybe this is -(l+ l-1 + ...)
            Ki[sum([nu-x for x in range(l)])+k_i[j] - l:sum([nu-x for x in range(l+1)]) ,shift+l] = L2
            Ki[sum([nu-x for x in range(j)])+k_i[l] - l:sum([nu-x for x in range(j+1)]) ,shift+l] = L1
    
    return(Ki)

# TODO: unfinished function
def comp_Kij_KKj(matrix, nu, m, n, r, k_ideal, i, j):
    field = matrix.base_ring()
    r_i = r[i]
    r_j = r[j]
    k_i = k_ideal[i]
    k_j = k_ideal[j]

    Kij_nrows = r_i*nu - r_i*(r_i-1)/2
    Kij_ncols = r_j*(r_j-1)/2
    Kij = zero_matrix(field, Kij_nrows, Kij_ncols)


    KKj_nrows = nu*(nu+1)/2 - r_i*nu + (r_i * (r_i-1)/2)
    KKj_ncols = r_j*(r_j-1)/2
    KKj = zero_matrix(field, KKj_nrows, KKj_ncols)
    col = 0

    for l in range(1,r_j):
        for m in range(l):
            S = comp_col_ij_h(matrix,nu,m,n,i,sum(r[:j])+l,k_j[m])
            new_row_order = koszul_order(nu, k_i)
            S_o = koszul_row_perm(S,new_row_order)
            Kij[:,col] = S_o[:Kij_nrows,:]
            KKj[:,col] = S_o[Kij_nrows:,:]
            col+=1
    return (Kij,KKj)
    
