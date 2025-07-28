from macaulay_matrix import *
from order import *
from code_ideal import *
from expanded_matrix import *

from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

import time


def extend(C, T, new_col_order):
    base = matrix([vector(list(-T*elt) + list(elt)) for elt in C])
    base = kosz_col_perm_inv(base, new_col_order)
    return base


class LinearStrand:
    """
    Linear strand
    """

    def __init__(self, nu, i, B_gen):
        """
        nu: number of variables
        i: index of B_gen
        B_gen: linear strand at range i
        """
        self.nu = nu
        self.i = i
        self.B_gen = B_gen
    

    def compute(self, r):
        base = self.B_gen
        index = self.i
        betti = []

        while True and index < r :
            betti.append(base.nrows())
            print("beta", base.nrows())
        
            start = time.time()
            B = ExpandedMatrix(self.nu, base)
            B.ef()

            T, TT, new_col_order = B.kosz_syz_mats()
            C = TT.right_kernel().matrix()
        
            end = time.time()
            print(f"[{round(end-start)}]loop")
        
            if C == 0 :
                break
        
            base = extend(C, T, new_col_order)
        
        return betti
    

class CodeStrand:
    """
    Linear strand of a code, syzygy distinguisher
    """


    def __init__(self, G):
        self.field = G.base_ring()
        self.G = G
        self.nu = G.nrows()


    def __repr__(self):
        pass


    def compute(self, r = None):
        if r == None :
            r = self.nu
        
        R = PolynomialRing(self.field, self.nu, 'x', order = 'lex')
        strand = [0]
    
        i2 = ideal2(self.G, R)
        strand.append(len(i2))
    
        M23 = macaulay23(i2, R)
        i3 = M23.right_kernel()

        B_gen = i3.matrix()
        syzygies = LinearStrand(self.nu, 2, B_gen)

        strand += syzygies.compute(r)

        return strand


def exemple_1():
    k = GF(11)
    nu = 5
    R = PolynomialRing(k, nu, 'x', order ="lex");
    x = R.gens()

    A = matrix (R, [
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

    #lin_b = lin_betti(B_gen, nu)
    #print("The linear strand is :")
    #print(lin_b)


def gaulay_11_6():
    G = matrix(GF(3),
    [[1, 0, 0, 0, 0, 0, 2, 0, 1, 2, 1],
    [0, 1, 0, 0, 0, 0, 1, 2, 2, 2, 1],
    [0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1],
    [0, 0, 0, 1, 0, 0, 1, 1, 0, 2, 2],
    [0, 0, 0, 0, 1, 0, 2, 1, 2, 2, 0],
    [0, 0, 0, 0, 0, 1, 0, 2, 1, 2, 2]])
    #code_betti(G, 6)
    golay = CodeStrand(G)
    dist = golay.compute(r=4)
    print(dist)

def hamming_7_4():
    G = matrix(GF(2),
    [[1, 0, 0, 0, 1, 1, 1],
    [0, 1, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1, 0, 0],
    [0, 0, 0, 1, 0, 0, 1]])
    #code_betti(G, 4)

    G = matrix(GF(2),
    [[1, 0, 0, 0, 1, 0, 1],
    [0, 1, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1, 1, 0],
    [0, 0, 0, 1, 0, 1, 1]])
    #code_betti(G, 4)

    
if __name__ == '__main__':
    gaulay_11_6()
