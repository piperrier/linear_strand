from macaulay_matrix import *
from order import *
from code_ideal import *
from expanded_matrix import *

from sage.matrix.constructor import matrix

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
        
            #start = time.time()
            B = ExpandedMatrix(self.nu, base)
            B.ef()

            T, TT, new_col_order = B.kosz_syz_mats()
            C = TT.right_kernel().matrix()
            #end = time.time()
        
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

        strand = [0]
        
        N = square_matrix(self.G).kernel().matrix()
        strand.append(N.nrows())

        M23 = macaulay23(N,self.nu)
        i3 = M23.right_kernel()

        B_gen = i3.matrix()
        syzygies = LinearStrand(self.nu, 2, B_gen)

        strand += syzygies.compute(r)

        return strand
