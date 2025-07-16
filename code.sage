# =====================
def carre_matrice(M):
    k = M.nrows()
    # lex
    L = [M.row(i).pairwise_product(M.row(j)) for i in range(0,k) for j in range(i,k)]
    return matrix(L)

# =====================
def cube_matrice(M):
    k = M.nrows()
    L = [M.row(i1).pairwise_product(M.row(i2)).pairwise_product(M.row(i3)) for i1 in range(0,k) for i2 in range(i1,k) for i3 in range(i2,k)]
    return matrix(L)    

# =====================
def ideal2(G,R):
    k = G.nrows()
    x = R.gens()
    
    N = carre_matrice(G).kernel()
    b12 = N.dimension()
    
    I = []
    if b12 >= 1:
        k2 = k*(k+1)/2
        M2 = [x[i]*x[j] for i in range(0,k) for j in range(i,k)]
        for i in range(b12):
            b = 0
            for j in range(k2):
                b+=N.basis()[i][j]*M2[j]
            I.append(b)
    
    return I

# =====================
def macaulay23(I, R):
    F = R.base_ring()

    x = R.gens()
    k = len(x)

    B3 = [x[i1]*x[i2]*x[i3] for i1 in range(k) for i2 in range(i1,k) for i3 in range(i2,k)]
    n = len(I)

    M = zero_matrix(R,len(B3),k*n)

    for i in range(n):
        cf = I[i].coefficients()
        mf = I[i].monomials()
        for j in range(k):
            for t in range(len(mf)):
                M[B3.index(mf[t]*x[j]),i*k+j] = cf[t]
    
    return M