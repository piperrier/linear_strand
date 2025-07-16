### ========================
# ideal row ordering
def kappa(nu, k_list):
    # assume k_list is ordered, from lowest to greatest
    kappa = []
    for j in range(len(k_list)) :
        sigma = [(i,k_list[j]) for i in range(0,k_list[j])] + [(k_list[j],i) for i in range(k_list[j], nu)]
        tau = [(k_list[i],k_list[j]) for i in range(0,j)]
        kappa += [elt for elt in sigma if elt not in tau]
    return kappa

def koszul_order(nu, k_list):
    order = UnorderedTuples(range(nu),2).list()

    kap = kappa(nu, k_list)
    new_order = kap + [elt for elt in order if elt not in kap]
    return new_order

# DONE
def koszul_row_perm(matrix, order):
    new_matrix = copy(matrix)
    for i, perm in enumerate(order) :
        shift = sum([i for i in range(nu, nu-perm[0],-1)])
        new_matrix[i] = matrix[shift+perm[1]-perm[0]]

    return new_matrix


### ============================
# module row ordering
def koszul_order_block(nu,m,n, k_ideal):
    order = UnorderedTuples(range(nu),2).list()
    block_order = [(i, elt) for i in range(0,m) for elt in order]
    l = []
    for i, k_i in enumerate(k_ideal):
        kappa_i = kappa(nu, k_i)
        l+= [x for x in zip([i]*len(kappa_i),kappa_i)]
        
    new_order = l + [elt for elt in block_order if elt not in l]
    return new_order

def block_kosz_row_perm(matrix, nu, order):
    new_matrix = copy(matrix)
    for i,(block,perm) in enumerate(order) :
        shift = sum([j for j in range(nu, nu-perm[0],-1)])
        new_matrix[i] = matrix[block * (nu*(nu+1)/2) + shift + perm[1]-perm[0]]
    return new_matrix

### ================================
# module col ordering
def lambdaa(r, k_ideal, nu):
    col = []
    for i, k_i in enumerate(k_ideal):
        if len(k_i)>1:
            for j in range(1,len(k_i)):
                for h in range(0,j):
                    col.append((sum(r[:i])+j)*nu + k_i[h])
    return col

def col_sol(r,k_ideal,nu):
    col = []
    for i,k_i in enumerate(k_ideal):
        if len(k_i)>1:
            for j in range(1,len(k_i)):
                for h in range(0,j):
                    col.append((sum(r[:i])+j,k_i[h]))
    return col

def col_order(r,k_ideal,nu,n):
    order = [i for i in range(nu*n)]
    lam = lambdaa(r,k_ideal,nu)
    new_order = [elt for elt in order if elt not in lam] + lam
    return new_order

def kosz_col_perm(matrix, order):
    new_matrix = copy(matrix)[:,order]
    return new_matrix

def kosz_col_perm_inv(matrix, order):
    new_matrix = copy(matrix)
    for i, perm in enumerate(order):
        new_matrix[:,perm] = matrix[:,i]
    return new_matrix

### =======================================
# unused (cause useless) and maybe false
def transvections(matrix,nu,m,n,r,k_ideal):
    # pour chaque bloque diag
    for i in range(len(r)):
        #print(i)
        r_i = r[i]
        k_i = k_ideal[i]
        # pour chaque pivot du bloque / chaque ligne de B / chaque bloque diag de Sij
        for j in range(r_i):
            pivot = k_i[j]
            #pour tous les coefs qui suivent le pivot
            for k in range(pivot+1,nu):
                ligne_B = sum(r[:i]) + j
                col_B = nu * i + k
                coeff = matrix[ligne_B,col_B]
                if coeff != 0 :
                #le supprimer pour toutes les cols qui sont multipliées par un "pivot"
                    for l in range(j,r_i):
                        shift = sum([nu-p for p in range(j)])
                        row_piv = shift+l
                        row_elim = l
                        print("shift", shift)
                        print("row piv", row_piv)
                        print("row elim", row_elim) 