from pylab import *
import numpy as np
import time


def jacobi(a,tol = 1.0e-9): # Jacobi method
    itera = 0.
    def maxElem(a): # Find largest off-diag. element a[k,l]
        n = len(a)
        aMax = 0.0
        for i in range(n-1):
            for j in range(i+1,n):
                if abs(a[i,j]) >= aMax:
                    aMax = abs(a[i,j])
                    k = i; l = j
        return aMax,k,l
 
    def rotate(a,p,k,l): # Rotate to make a[k,l] = 0
        n = len(a)
        aDiff = a[l,l] - a[k,k]
        #print 'Differense all-akk=',aDiff
        if abs(a[k,l]) < abs(aDiff)*1.0e-36:
            t = a[k,l]/aDiff
            #print 'abs(a[k,l]) < abs(aDiff)', t
        else:
            phi = aDiff/(2.0*a[k,l])
            t = 1.0/(abs(phi) + sqrt(phi**2 + 1.0))
            #print 't: ',t
            if phi < 0.0: t = -t
        c = 1.0/sqrt(t**2 + 1.0)
        s = t*c
        tau = s/(1.0 + c)
        temp = a[k,l]
        #print 'k,l ', k,l
        #print 'values of point akl,akk,all\n',a[k,l],a[k,k],a[l,l]
        a[k,l] = 0.0
        a[k,k] = a[k,k] - t*temp
        a[l,l] = a[l,l] + t*temp

        
        for i in range(k):      # Case of i < k
            temp = a[i,k]
            a[i,k] = temp - s*(a[i,l] + tau*temp)
            a[i,l] = a[i,l] + s*(temp - tau*a[i,l])
        for i in range(k+1,l):  # Case of k < i < l
            temp = a[k,i]
            a[k,i] = temp - s*(a[i,l] + tau*a[k,i])
            a[i,l] = a[i,l] + s*(temp - tau*a[i,l])
        for i in range(l+1,n):  # Case of i > l
            temp = a[k,i]
            a[k,i] = temp - s*(a[l,i] + tau*temp)
            a[l,i] = a[l,i] + s*(temp - tau*a[l,i])
        for i in range(n):      # Update transformation matrix
            temp = p[i,k]
            p[i,k] = temp - s*(p[i,l] + tau*p[i,k])
            p[i,l] = p[i,l] + s*(temp - tau*p[i,l])
        #print 'Matrix after rot\n',a
    n = len(a)
    maxRot = 5*(n**2)       # Set limit on number of rotations
    p = identity(n)*1.0 
    for i in range(maxRot):
        itera += 1 
        aMax,k,l = maxElem(a)
        if aMax < tol: return diagonal(a),p,itera
        rotate(a,p,k,l)

def A(n):
    rho_min = 0.
    rho_max = 5.
    
    #n = 50
    A = zeros(shape=(n-1, n-1))
    
    rho = zeros(n+1)
    rho[0] = rho_min
    V = zeros(n-1)
    d = zeros(n-1)
    
    h = (rho_max - rho_min)/n
    e = -1./(h**2);
    
    for i in range(n+1):
        rho[i] = rho[0] + i*h
    
    for i in range(n-1):
    
        V[i] = rho[i+1]**2 # Non-interaction case
        d[i] = 2./(h**2) + V[i]

    #Initialize A
    A[n-2,n-2] = d[n-2]
    for i in range(n-2):
        A[i,i] = d[i]
        A[i,i+1] = e
        A[i+1,i] = e

    return A,h,rho
#print 'Start matrix\n',A


print '_n__________________eigenvalues_________________________time[s]______iterations______'
n = linspace(50,175,6)
for n in n:
    n = int(n)
    a,h,rho = A(n)
    t0 = time.clock()
    eigenval,eigenvec,itera = jacobi(a)
    t1 = time.clock() - t0
    permute = eigenval.argsort()
    eigenval = eigenval[permute]
    print '{:<6}{:<46}{:<20}{:<10}'.format(n,eigenval[0:3],t1,itera)
    

t0 = time.clock()
a,h,rho = A(175)
eigenval,eigenvec,itera = jacobi(a)
t1 = time.clock() - t0

permute = eigenval.argsort()
eigenval = eigenval[permute]
eigenvec = eigenvec[:,permute]


eigenvec1 = eigenvec[:,0]
eigenvec2 = eigenvec[:,1]
eigenvec3 = eigenvec[:,2]

norm1 = 0
norm2 = 0
norm3 = 0

u_square_norm1 = zeros(len(eigenvec)-1)
u_square_norm2 = zeros(len(eigenvec)-1)
u_square_norm3 = zeros(len(eigenvec)-1)


norm1 = sum(eigenvec1**2)*h
norm2 = sum(eigenvec2**2)*h
norm3 = sum(eigenvec3**2)*h

u_square_norm1 = eigenvec1**2 / norm1
u_square_norm2 = eigenvec2**2 / norm2
u_square_norm3 = eigenvec3**2 / norm3


figure()
plot(rho[0:-2],u_square_norm1,'-r',
     rho[0:-2],u_square_norm2,'-g',
     rho[0:-2],u_square_norm3,'-b',)
grid(True)
title('Wavefunction for n={}'.format(n))
xlabel(r'$\rho$')
ylabel(r'$|u(\rho)|^2$')
savefig('non-interaction_case.png')
show()
