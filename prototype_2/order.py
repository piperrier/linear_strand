from sage.combinat.tuple import UnorderedTuples
from sage.combinat.permutation import Permutation


#################################
### Rows order
#################################


def kappa(nu, k_list):
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


def koszul_order_block(nu, m, n, k_ideal):
    order = UnorderedTuples(range(nu),2).list()
    block_order = [(i, elt) for i in range(0,m) for elt in order]
    l = []
    
    for block_i, k_i in enumerate(k_ideal):
        kappa_i = kappa(nu, k_i)
        l += [x for x in zip([block_i]*len(kappa_i),kappa_i)]
        
    new_order = l + [elt for elt in block_order if elt not in l]
    
    return new_order


def block_kosz_row_perm(matrix, nu, order):
    ord_to_int = lambda block, perm : (block+1) * (nu * (nu + 1)//2) - ((nu - perm[0]) * (nu - perm[0] + 1) // 2) + perm[1] - perm[0]
    new_order = [ord_to_int(block, perm) for (block, perm) in order]
    new_matrix = matrix[new_order,:]
    
    return new_matrix


#################################
### Columns order
#################################


def lambdaa(r, k_ideal, nu):
    col = []

    for block_i, k_i in enumerate(k_ideal):
        if len(k_i)>1:
            for j in range(1,len(k_i)):
                for h in range(0,j):
                    col.append((sum(r[:block_i])+j)*nu + k_i[h])
    
    return col


def col_order(r, k_ideal, nu, n):
    order = [i for i in range(nu*n)]
    lam = lambdaa(r,k_ideal,nu)
    new_order = [elt for elt in order if elt not in lam] + lam
    
    return new_order


def kosz_col_perm(matrix, order):
    new_matrix = matrix[:,order]
    return new_matrix


def kosz_col_perm_inv(matrix, order):
    order = [x + 1 for x in order]
    perm = [x - 1 for x in list(Permutation(order).inverse())]        
    new_matrix = matrix[:,perm]
    
    return new_matrix