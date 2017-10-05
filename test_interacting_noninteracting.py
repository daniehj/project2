from interacting import A
from pylab import identity,array


n = 5
rho_max = 5
omega =1
A,h,rho = A(n,omega,rho_max)
#Check if sum of the diagonal in array A
#is the same as the sum of h/2 + V
def get_diag(A):
    I = identity(len(A))
    diag = A*I
    
    h =  sum(A*I)
    diag =  sum(h)
    return diag

def get_d(rho):

    V = omega**2*rho[1:]**2 + 1/rho[1:]
    d = 2./(h**2) + V
    d = d[:-1]
    return sum(d)

def test_diagonal():
    assert get_diag(A) == get_d(rho)
