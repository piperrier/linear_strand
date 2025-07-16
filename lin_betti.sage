load("utils.sage")
load("matrix.sage")
load("combinatorics.sage")
load("code.sage")

import time

def kosz_syz_mats(B, nu):

    #pre_processing
    n,m,r,k_ideal,t,t_details = pre_processing(B,nu)

    S = matrix_syzygetic(B,nu,m,n)

    #row order
    new_row_order = koszul_order_block(nu,m,n,k_ideal)
    S_row_perm = block_kosz_row_perm(S,nu,new_row_order)

    #col order
    new_col_order = col_order(r,k_ideal,nu,n)
    S_col_perm = kosz_col_perm(S_row_perm,new_col_order)
    S_col_perm = t_subdivide(S_col_perm,t)

    #ideal elim
    S_ie = block_ideal_elim(S_col_perm,r,k_ideal, t_details)
    
    #block_forward
    S_bf = block_forward_substitution(S_ie)

    # T,TT
    T,TT = block_T_TT(S_bf)
    return (T,TT,new_col_order)

def extend(C,T, new_col_order):
    # TODO new col order
    base = matrix([vector(list(-T*elt) + list(elt)) for elt in C])
    base = kosz_col_perm_inv(base, new_col_order)
    return base

def lin_betti(B_gen, nu):
    # TODO
    
    betti = [B_gen.nrows()]
    print("beta_23", B_gen.nrows())

    start = time.time()
    
    B = matrix_rru(B_gen)
    T,TT, new_col_order = kosz_syz_mats(B, nu)
    #C = TT.right_kernel().basis()
    C = TT.right_kernel().matrix()

    end = time.time()
    print(f"[{end-start}]lin_betti")


    #while C != [] :
    while C != 0 :
        start = time.time()

        B = extend(C,T, new_col_order)
        betti.append(B.nrows())
        print("beta", B.nrows())

        B = matrix_rru(B)
        T,TT, new_col_order = kosz_syz_mats(B,nu)
        #C = TT.right_kernel().basis()
        C = TT.right_kernel().matrix()
        
        end = time.time()
        print(f"[{end-start}]loop")

    return betti

def code_betti(G_gen, nu):
    k = G_gen.base_ring()
    
    R = PolynomialRing(k, nu, 'x', order ="lex");
    x = R.gens()

    lin_b = [0]

    start = time.time()
    i2 = ideal2(G_gen,R)
    end = time.time()
    lin_b.append(len(i2))
    print(f"[{end-start}]beta_12:{len(i2)}")


    M23 = macaulay23(i2, R)

    start = time.time()
    i3 = M23.right_kernel()
    end = time.time()
    print(f"[{end-start}]M23")

    B_gen = i3.matrix()
    lin_b += lin_betti(B_gen,nu)

    print("The linear strand is :")
    print(lin_b)

def exemple_1():
    k = GF(11)
    nu = 5
    R = PolynomialRing(k, nu, 'x', order ="lex");
    x = R.gens()

    A = Matrix (R, [
    [x[2],x[3],x[4],0,0,0,0,0],
    [0,0,0,x[1],x[3],x[4],0,0],
    [x[1]+9*x[2],7*x[1],2*x[1]+10*x[4],5*x[1]+x[2],2*x[2],9*x[2]+6*x[4],0,0],
    [10*x[0],0,0,10*x[0],0,0,x[3],x[4]],
    [10*x[0]+6*x[2],3*x[0],9*x[0]+3*x[4],6*x[0]+8*x[2],0,0,x[2],10*x[2] + 3*x[4]],
    [2*x[0] + 5*x[1], 0, 0,10*x[0] + 3*x[1],8*x[0],2*x[0]+8*x[4],9*x[1],x[1]+5*x[4]]
    ])
    
    print(A)
    print("\n")

    B_gen = construct_B_from_A(A)
    print(B_gen)
    print("\n")

    lin_b = lin_betti(B_gen, nu)
    print("The linear strand is :")
    print(lin_b)

def gaulay_11_8():
    G = matrix(GF(3),
    [[1, 0, 0, 0, 0, 0, 2, 0, 1, 2, 1],
    [0, 1, 0, 0, 0, 0, 1, 2, 2, 2, 1],
    [0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1],
    [0, 0, 0, 1, 0, 0, 1, 1, 0, 2, 2],
    [0, 0, 0, 0, 1, 0, 2, 1, 2, 2, 0],
    [0, 0, 0, 0, 0, 1, 0, 2, 1, 2, 2]])

    code_betti(G,6)
    
def hamming_7_4():
    G = matrix(GF(2),
    [[1, 0, 0, 0, 1, 1, 1],
    [0, 1, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1, 0, 0],
    [0, 0, 0, 1, 0, 0, 1]])
    code_betti(G,4)

    G = matrix(GF(2),
    [[1, 0, 0, 0, 1, 0, 1],
    [0, 1, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1, 1, 0],
    [0, 0, 0, 1, 0, 1, 1]])
    code_betti(G,4)

    
if __name__ == '__main__' and '__file__' in globals():
    
    gaulay_11_8()
