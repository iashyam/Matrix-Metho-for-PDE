import matplotlib.pyplot as plt 
import numpy as np
from scipy.linalg import det
import cmath
import math

#<== Constants ==>
h = 1 #plank constant
m = 1 #mass of the object
v0 = 1 #barrier height
a = 1
i = 1j #defining iota 
E = np.linspace(0, 2, 200)


Rs = []  # list of reflecing and trasmattting constant
Ts = []
total = []  # list of r+t
T_th = []

#List barriers as tuple if starting point, endpoint 
# we are doing this because barrier height is same for all
barriers = [(0,a), (2*a, 3*a), (5*a, 7*a)]


def scatterinMatrix(k1, k2):
    d = np.zeros([2, 2], dtype=complex)
    d[0][0] = (1+k1/k2)
    d[0][1] = (1-k1/k2)
    d[1][1] = (1+k1/k2)
    d [1][0] = (1-k1/k2)

    return (1/2)*d

def progogationMatrix(a, k):

    return np.array([[cmath.exp(i*k*a), 0], [0, cmath.exp(-i*k*a)]], dtype=complex)


# def matrix_method(k1,k2,k3, n):
#     a,b = barriers[0]
#     c,d = barriers[1]
#     d21 = scatterinMatrix(k1, k2)
#     d32 = scatterinMatrix(k2, k3)

#     p1 = progogationMatrix(a, k1)
#     p2 = progogationMatrix(b, k2)
#     p3 = progogationMatrix(c, k1)
#     p4 = progogationMatrix(d, k2)


#     m1 = np.matmul(d21, p3)
#     m2 = np.matmul(d32, p4)
#     m = np.matmul(m2,m1)

#     M1 = np.matmul(d21, p1)h
#     M2 = np.matmul(d32, p2)
#     M =  np.matmul(M2, M1)

#     return np.matmul(m, M)


'''
The total effect of two metrices are basically 
(d3<-2)(P(a,k)) times this metrix for one potential box
'''
def matrix_method(k1,k2,k3, n):
    a,b = barriers[n-1]

    d21 = scatterinMatrix(k1, k2)
    d32 = scatterinMatrix(k2, k3)

    p1 = progogationMatrix(a, k1)
    p2 = progogationMatrix(b, k2)

    M1 = np.matmul(d21, p1)
    M2 = np.matmul(d32, p2)
    M =  np.matmul(M2, M1)

    if n == 1:
        return M
    else:
        return np.matmul(M,matrix_method(k1,k2,k3,n-1))

def theory(e):
    return 16*(e/v0)*math.exp((-2*a*math.sqrt(2*m*v0)/h**2)**2)
    
for e in E:
    k1 = cmath.sqrt(2*m*e/(h**2))
    k2 = cmath.sqrt((2*m*(e-v0))/(h**2))
    k3 = k1

    try:
        M = matrix_method(k1,k2,k3,2)
        # print(M)

        r = -(M[1][0]/M[1][1])
        t = det(M)/M[1][1]

        R = r*np.conjugate(r)
        T = t*np.conjugate(t)

        Rs.append(R)
        Ts.append(T)

        T_th.append(theory(e))
       
        total.append(T+R)
    except ZeroDivisionError:
         E = np.setdiff1d(E, np.array(e))

# plt.plot(T_th,Ts , label='Theory')
plt.plot(E, Rs, label="Reflection")
plt.plot(E, Ts, label="Transmission")
plt.plot(E, total)
plt.ylabel('Energy')
plt.title('Energy vs Transmission and Reflective Cofficients')
plt.xlabel('R, T')
plt.legend()
plt.show()