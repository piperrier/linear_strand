from sage.matrix.constructor import matrix
from sage.matrix.special import zero_matrix
from sage.functions.other import binomial


def square_matrix(M):
    k = M.nrows()
    L = [M.row(i).pairwise_product(M.row(j)) for i in range(0,k) for j in range(i,k)]
    return matrix(L)


def index_to_mon_deg2(index, k):
    '''
    associe à un indice le monôme de degré 2 à k variables associé dans l'ordre lex
    '''
    mon = []
    x=0
    i = -1
    while x <= index :
        i+=1
        x+=k-i

    mon.append(i)
    x -= (k-i)
    mon.append(index-x+i)
    return mon


def mon_deg3_to_index(mon, k):
    '''
    associe à un monôme de degré 3 l'indice associé dans la liste
    '''
    index = 0

    i0 = mon[0]
    index += binomial(k+2,3) - binomial(k+1-(i0-1),3)

    i1 = mon[1]
    index += (k-i0)*(k-i0+1)//2 - (k-i1+1)*(k-i1)//2

    i2 = mon[2]
    index += i2 - i1

    return index


def macaulay23(N, k):
    b12 = N.nrows()

    M = zero_matrix(N.base_ring(), binomial(k+2,3)  ,k*b12)

    for i in range(b12):
        d = N.row(i).dict()
        for j in range(k):
        #À i est associé un monôme de degré 2 x_l*x_p, on lui associe x_l*x_p*x_j et on construit la colonne associée
            for index_mon2,coef in d.items():
                mon = index_to_mon_deg2(index_mon2,k)
                mon.append(j)
                mon.sort()

                M[mon_deg3_to_index(mon,k),i*k+j] = coef

    return M
