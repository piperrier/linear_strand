K = GF(11)
nu = 5
R = PolynomialRing(K,'x',nu)
x = R.gens()


#convention pour B e_0,...,e_m pour l'ordre des colonnes
def construct_A_from_B(B):
    m = B.ncols()//nu
    n = B.nrows()
    A = matrix(R,m,n)
    for i in range(n):
        for j in range(m):
            for p in range(nu):
                A[j,i] += B[i,j*nu+p]*x[p]
    return A




def construct_B_from_A(A):
    m = A.nrows()
    n = A.ncols()
    B = matrix(K,n,nu*m)
    for i in range(n):
        for j in range(m):
            for p in range(nu):
                B[i,j*nu +p] += A[j,i].monomial_coefficient(x[p])
    return B


#on suppose B RLU

def construct_k_from_B(B):
    n=B.nrows()
    m=B.ncols()//nu
    k=[]
    r=[]
    for i in range(m):
        block = B[sum(r):,i*nu:(i+1)*nu]
        row_index=0
        k_i=[]
        while row_index <n-sum(r) and block.row(row_index)!=0:
            k_i.append(block.row(row_index).trailing_support())
            row_index +=1
        r.append(row_index)
        k.append(k_i)
    return k

    


def construct_K(B,k):
    m = B.ncols()//nu
    n = B.nrows()

    #Calcul des dimensions de T
    ncol = 0
    for i in range(m):
        ncol += len(k[i])*(len(k[i])-1)//2
    nrow = n*nu - ncol

    #calcul de sigma :
    s = 0
    sigma = [0]
    for i in range(m):
        s += len(k[i])
        sigma.append(s)
    #calcul du début des lignes et colonnes des blocs
    l = 0
    c=0
    colonne_debut_bloc = [0]
    ligne_debut_bloc = [0]
    for i in range(m):
        l += len(k[i])*nu - len(k[i])*(len(k[i])-1)//2
        ligne_debut_bloc.append(l)
        c += len(k[i])*(len(k[i])-1)//2
        colonne_debut_bloc.append(c)


    
    #On construit T0 qui correspond aux colonnes solutions ordonnés comme dans S puis on fera la permutation de lignes
    T0 = Matrix(K,m*nu*(nu+1)//2,ncol)
    #on construit les K_i
    colonne = 0
    for i in range(m):
        for q in range(len(k[i])):
            for l in range(q):
                #on construit la colonne x_k[i][l]*L_i,sigma[i]+q
                
                #on construit L_i,sigma[i]+q au niveau du bloc k_l
                for p in range(k[i][l],nu):
                    T0[i*nu*(nu+1)//2 + k[i][l]*nu - k[i][l]*(k[i][l]-1)//2 + p-k[i][l],colonne] = B[sigma[i]+q,i*nu+p]

                #on construit -L_l au niveau du bloc k_q
                for p in range(k[i][l]+1,k[i][q]):
                    T0[i*nu*(nu+1)//2 + p*nu-p*(p-1)//2 + k[i][q]-p,colonne] = -1*B[sigma[i]+l,i*nu+p]
                # skip step : p = k[i][q]
                for p in range(k[i][q]+1,nu):
                    T0[i*nu*(nu+1)//2 +   k[i][q]*nu-k[i][q]*(k[i][q]-1)//2 + p-k[i][q],colonne] = -1*B[sigma[i]+l,i*nu+p]
                    

                colonne += 1
    
    #on construit les K_ij:
    for i in range(1,m):
        colonne = 0
        for j in range(i):
            for q in range(len(k[j])):
                for l in range(q):
                    #on construit la colonne x_k[j][l]*L_i,(sigma[j]+q)
                    for p in range(k[j][l]):
                        T0[i*nu*(nu+1)//2 + p*nu-p*(p-1)//2 + k[j][l]-p,colonne] = B[sigma[j]+q,i*nu+p]
                    for p in range(k[j][l],nu):
                        T0[i*nu*(nu+1)//2 + k[j][l]*nu-k[j][l]*(k[j][l]-1)//2 + p-k[j][l],colonne] = B[sigma[j]+q,i*nu+p]
                    colonne +=1

    T = matrix(K,nrow,ncol)
    Tbar = matrix(K,m*nu*(nu+1)//2 - nrow,ncol)


    #on réarrange les lignes de T0 pour avoir T et Tbar
    ligne = 0
    for i in range(m):
        for q in range(len(k[i])):
            for p in range(nu):
                if p not in k[i][:q]:
                    #on construit la ligne x_p*x_k[i][q]
                    mk = min(p,k[i][q])
                    T[ligne] = T0[i*nu*(nu+1)//2 + mk*nu-mk*(mk-1)//2 + max(p,k[i][q])-mk]
                    ligne +=1
    ligne=0
    for i in range(m):
        for q in range(nu):
            if q not in k[i]:
                for p in range(q,nu):
                    if p not in k[i]:
                        Tbar[ligne] = T0[i*nu*(nu+1)//2 + q*nu-q*(q-1)//2 + p-q]
                        ligne +=1



    return(T,Tbar)

#renvoie H_ij plus un autre morceau qui correspond aux blocs de ligne de Hbar_j correspondant à i
def construct_H(B,k,i,j):
    m = B.ncols()//nu
    n = B.nrows()

    #calcul de sigma :
    s = 0
    sigma = [0]
    for l in range(m):
        s += len(k[l])
        sigma.append(s)
    
    #H0 est S_ij sans les colonnes solutions et dans l'ordre des lignes original
    H0 = Matrix(K,nu*(nu+1)//2,len(k[j])*nu-len(k[j])*(len(k[j])-1)//2)
    colonne =0
    for q in range(len(k[j])):
        for p in range(nu):
            if p not in k[j][:q]:
                #on construit la colonne x_p*L_i,sigma[j]+q
                for l in range(nu):
                    mk = min(l,p)
                    H0[mk*nu - mk*(mk-1)//2 + max(l,p)-mk,colonne] = B[sigma[j]+q,i*nu+l]
                colonne += 1
    
    #on réordonne les lignes de H0 selon celle de l'idéal I_i
    H = Matrix(K,len(k[i])*nu-len(k[i])*(len(k[i])-1)//2,len(k[j])*nu-len(k[j])*(len(k[j])-1)//2)
    Hbar = Matrix(K,nu*(nu+1)//2 - len(k[i])*nu+len(k[i])*(len(k[i])-1)//2,len(k[j])*nu-len(k[j])*(len(k[j])-1)//2)
    
    ligne = 0
    for q in range(len(k[i])):
        for p in range(nu):
            if p not in k[i][:q]:
                #on construit la ligne x_p*x_k[i][q]
                mk = min(k[i][q],p)
                H[ligne] = H0[mk*nu-mk*(mk-1)//2 + max(k[i][q],p)-mk]
                ligne +=1
    ligne = 0
    for q in range(nu):
        if q not in k[i]:
            for p in range(q,nu):
                if p not in k[i]:
                    Hbar[ligne] = H0[q*nu-q*(q-1)//2 +p-q]
                    ligne += 1
    return (H,Hbar)


def construct_T(B,k):
    (T,Tbar) = construct_K(B,k)
    n=B.nrows()
    m=B.ncols()//nu

    debut_bloc_Hbar = [0]
    nrow_Hbar = 0
    for i in range(m):
        nrow_Hbar += nu*(nu+1)//2 - len(k[i])*nu+len(k[i])*(len(k[i])-1)//2
        debut_bloc_Hbar.append(nrow_Hbar)

    
    c = 0
    debut_bloc = [0]
    for i in range(m):
        c += len(k[i])*nu-len(k[i])*(len(k[i])-1)//2
        debut_bloc.append(c)
    
    for j in range(m):
        Hbar_j = Matrix(K,debut_bloc_Hbar[j],len(k[j])*nu-len(k[j])*(len(k[j])-1)//2)
        (_,Hbar) = construct_H(B,k,j,j)
        Hbar_j = Hbar_j.stack(Hbar)

        for i in range(j+1,m):
            (H,Hbar) = construct_H(B,k,i,j)
            T[debut_bloc[i]:debut_bloc[i+1]] = T[debut_bloc[i]:debut_bloc[i+1]] - H*T[debut_bloc[j]:debut_bloc[j+1]]
            Hbar_j = Hbar_j.stack(Hbar)

        #on simplifie Hbar_j
        Tbar = Tbar - Hbar_j*T[debut_bloc[j]:debut_bloc[j+1]]

    return(T,Tbar)


def extend(C,T,k):
    #on étend les vecteurs de C à l'aide de T puis on les remet dans le bon ordre
    base0 = matrix([vector(list(-T*elt) + list(elt)) for elt in C])
    base = matrix(K,base0.nrows(),base0.ncols())
    m = len(k)

    sigma = [0]
    c=0
    for i in range(m):
        c+=len(k[i])
        sigma.append(c)
    
    #on range les colonnes de base0 dans le bon ordre dans base
    colonne =0
    for i in range(m):
        for l in range(len(k[i])):
            for p in range(nu):
                if p not in k[i][:l]:
                    #on construit la colonne x_p*L_i,sigma[i]+l
                    base[:,sigma[i]*nu + l*nu +p] = base0[:,colonne]
                    colonne +=1
    for i in range(m):
        for j in range(1,len(k[i])):
            for l in range(j):
                #on construit la colonne correspondant à x_k[i][l]*L_i,(sigma[i]+j)
                base[:,sigma[i]*nu + j*nu + k[i][l]] = base0[:,colonne]
                #print(i,j,k[i][l],sigma[i]*nu + j*nu + k[i][l])
                colonne +=1

    return base.echelon_form()


def lin_betti(B_gen):
    betti = [B_gen.nrows()]
    print("beta_23", B_gen.nrows())

    B = B_gen.echelon_form()
    k = construct_k_from_B(B)
    (T,Tbar) = construct_T(B,k)
    C = Tbar.right_kernel().matrix()
    print(f"C\n{C}")

    while C.nrows()>0:
        B = extend(C,T,k)
        print(f"B_extend\n{B}")

        betti.append(B.nrows())
        print(B.nrows())

        k = construct_k_from_B(B)
        (T,Tbar) = construct_T(B,k)
        C = Tbar.right_kernel().matrix()
        print(f"C\n{C}")
    return betti



A = Matrix (R, [
[x[2],x[3],x[4],0,0,0,0,0],
[0,0,0,x[1],x[3],x[4],0,0],
[x[1]+9*x[2],7*x[1],2*x[1]+10*x[4],5*x[1]+x[2],2*x[2],9*x[2]+6*x[4],0,0],
[10*x[0],0,0,10*x[0],0,0,x[3],x[4]],
[10*x[0]+6*x[2],3*x[0],9*x[0]+3*x[4],6*x[0]+8*x[2],0,0,x[2],10*x[2] + 3*x[4]],
[2*x[0] + 5*x[1], 0, 0,10*x[0] + 3*x[1],8*x[0],2*x[0]+8*x[4],9*x[1],x[1]+5*x[4]]
])

print(A)


#A = Matrix (R, [[x[0],x[1],x[2],0,0,0,0,0], [0,0,0,x[0],x[1],x[3],0,0],[-1*x[0]+2*x[3],-4*x[3],-2*x[2]+x[3],-5*x[0],2*x[2],x[2]+5*x[3],0,0],[0,0,-1*x[4],0,0,-1*x[4],x[0],x[1]],
#[3*x[0]-2*x[4],3*x[4],-5*x[2]+2*x[4],-3*x[0]+2*x[4],-3*x[4],3*x[3]-x[4],5*x[0]+x[3],-2*x[3]],[0,0,5*x[3]+2*x[4],-3*x[0]+2*x[4], -3*x[4],3*x[3]-x[4], 5*x[0]+x[3],-2*x[3]]])
#print(A)

B = construct_B_from_A(A)

lin_betti(B)


