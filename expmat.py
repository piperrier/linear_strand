from mac_matrix import *
from order import *

class ExpandedMatrix:
    """
    This class represent the expanded polynomial matrix
    """

    
    def __init__(self, nu, matrix):
        self.nu = nu
        self.n = matrix.nrows()
        self.m = matrix.ncols()//nu
        self.matrix = matrix


    def ef(self):
        #self.matrix.echelonize()       
        self.matrix = self.matrix.echelon_form()
    

    def pre_processing(self):
        r = []
        k = []
        for i in range(0,self.m):
            block = self.matrix[sum(r):,i*self.nu:(i+1)*self.nu]
            row_index = 0
            k_i = []
            while row_index < self.n-sum(r) and block.row(row_index) != 0 :
                k_i.append(block.row(row_index).trailing_support())
                row_index += 1
            r.append(row_index)
            k.append(k_i)
        t = [self.nu*ri-(ri*(ri-1))//2 for ri in r]
        spacing = [[self.nu-x for x in range(ri)] for ri in r]
        
        self.r = r
        self.k = k
        self.t = t
        self.spacing = spacing
        self.n_sol_col = [(ri*(ri-1))//2 for ri in self.r]


    def __repr__(self):
        return "Expanded matrix in %s :\n%s" % (self.matrix.base_ring(), self.matrix)


    def settings(self):
        print(f"nu={self.nu}")
        print(f"n={self.n}")
        print(f"m={self.m}")
        print(f"nu*(nu+1)/2={self.nu*(self.nu+1)//2}")
        print(f"r={self.r}")
        print(f"k={self.k}")
        print(f"t={self.t}")
        print(f"spacing={self.spacing}")
        print(f"n_sol_col={self.n_sol_col}")


    def kosz_syz_mats(self):
        self.pre_processing()

        S = matrix_syzygetic(self.matrix, self.nu, self.m, self.n)

        #row order
        new_row_order = koszul_order_block(self.nu, self.m, self.n, self.k)
        S_row_perm = block_kosz_row_perm(S, self.nu, new_row_order)

        #col order
        new_col_order = col_order(self.r, self.k, self.nu, self.n)
        S_col_perm = kosz_col_perm(S_row_perm, new_col_order)
        S_col_perm = t_subdivide(S_col_perm, self.t)

        #ideal elim
        S_ie = block_ideal_elim(S_col_perm, self.r, self.k, self.spacing)
    
        #block_forward
        S_bf = block_forward_substitution(S_ie)

        # T,TT
        T, TT = block_T_TT(S_bf)
        return (T, TT, new_col_order)