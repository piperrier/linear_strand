class ExpMatrix:
    def __init__(self, nu, matrix):
        self.nu = nu
        self.n = matrix.nrows()
        self.m = matrix.ncols()/nu
        self.matrix = matrix
    
    def ef(self):
        m_rru = self.matrix.echelon_form()
        self.matrix = m_rru
    
    def pre_processing(self):
        r = []
        k = []
        for i in range(0,self.m):
            block = self.matrix[sum(r):,i*self.nu:(i+1)*self.nu]
            row_index = 0
            k_i = []
            while row_index< self.n-sum(r) and block.row(row_index)!=0 :
                k_i.append(block.row(row_index).trailing_support())
                row_index +=1
            r.append(row_index)
            k.append(k_i)
        t = [self.nu*ri-(ri*(ri-1))/2 for ri in r]
        spacing=[[self.nu-x for x in range(ri)] for ri in r]
        
        self.r = r
        self.k = k
        self.t = t
        self.spacing = spacing
        self.n_sol_col = [(ri*(ri-1))/2 for ri in self.r]

    def __repr__(self):
        return "Expanded matrix in %s :\n%s" % (self.matrix.base_ring(), self.matrix)

    def settings(self):
        print(f"nu={self.nu}")
        print(f"n={self.n}")
        print(f"m={self.m}")
        print(f"nu*(nu+1)/2={self.nu*(self.nu+1)/2}")
        print(f"r={self.r}")
        print(f"k={self.k}")
        print(f"t={[self.nu*ri-(ri*(ri-1))/2 for ri in self.r]}")
        print(f"spacing={self.spacing}")
        print(f"n_sol_col={self.n_sol_col}")
    
    #def kosz_syz_mats(self):