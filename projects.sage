from linear_strand import *


def gaulay_11_6():
    G = matrix(GF(3),
    [[1, 0, 0, 0, 0, 0, 2, 0, 1, 2, 1],
    [0, 1, 0, 0, 0, 0, 1, 2, 2, 2, 1],
    [0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1],
    [0, 0, 0, 1, 0, 0, 1, 1, 0, 2, 2],
    [0, 0, 0, 0, 1, 0, 2, 1, 2, 2, 0],
    [0, 0, 0, 0, 0, 1, 0, 2, 1, 2, 2]])
    golay = CodeStrand(G)
    dist = golay.compute(r=4)
    print(dist)

def hamming_7_4():
    G = matrix(GF(2),
    [[1, 0, 0, 0, 1, 1, 1],
    [0, 1, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1, 0, 0],
    [0, 0, 0, 1, 0, 0, 1]])
    ham1 = CodeStrand(G)
    dist1 = ham1.compute(r=4)
    print(dist1)

    G = matrix(GF(2),
    [[1, 0, 0, 0, 1, 0, 1],
    [0, 1, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1, 1, 0],
    [0, 0, 0, 1, 0, 1, 1]])
    ham2 = CodeStrand(G)
    dist2 = ham2.compute(r=4)
    print(dist2)

    
if __name__ == '__main__':
    gaulay_11_6()