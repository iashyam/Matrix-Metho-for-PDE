import matplotlib.pyplot as plt 
import numpy as np
from scipy.linalg import det
import cmath
import math

#<== Constants ==>
h = 1 #plank constant
m = 1 #mass of the object
v0 = 0.3 #barrier height
# a = 1
i = 1j #defining iota 
E = np.linspace(0, 0.8, 400)


Rs = []  # list of reflecing and trasmattting constant
Ts = []
total = []  # list of r+t
T_th = []

#List barriers as tuple if starting point, endpoint 
# we are doing this because barrier height is same for all
barriers = [(4,8),(13,17)]


# def get_barriers(number, width, distance):
#     barriers = [(0, width)]
#     for n in range(1, number):
#         barriers.append((0+n*distance, width+n*distance))
#     return barriers

# v0 = eval(input('Enter the barrier height: '))
# number = eval(input("Enter number of potential: "))
# width = eval(input("Enter the width of one barrier: "))
# distance = eval(input('Enter the distance between two barriers: '))
# barriers = get_barriers(number, width, distance)



def scatterinMatrix(k1, k2):
    d = np.zeros([2, 2], dtype=complex)
    d[0][0] = (1+k1/k2)
    d[0][1] = (1-k1/k2)
    d[1][1] = (1+k1/k2)
    d [1][0] = (1-k1/k2)

    return (1/2)*d

def progogationMatrix(a, k):

    return np.array([[cmath.exp(i*k*a),0], [0,cmath.exp(-i*k*a)]], dtype=complex)


'''
The total effect of two metrices are
(d3<-2)(P(a,k)) times this metrix for one potential box
'''
def matrix_method(k1,k2,k3, n):
    a,b = barriers[n-1]

    d21 = scatterinMatrix(k1, k2)
    d32 = scatterinMatrix(k2, k3)

    p1 = progogationMatrix(a, k1)
    p2 = progogationMatrix(b, k2)

    M1 = np.matmul(p1, d21)
    M2 = np.matmul(p2, d32)
    M =  np.matmul(M1, M2)

    if n == 1:
        return M
    else:
        return np.matmul(M,matrix_method(k1,k2,k3,n-1))

# def func(e):
#     k1 = cmath.sqrt(2*m*e/(h**2))
#     k2 = cmath.sqrt((2*m*(e-v0))/(h**2))
#     k3 = k1

#     d21 = scatterinMatrix(k1, k2)
#     d32 = scatterinMatrix(k2, k3)

#     print(d21, d32)

# func(0.3)

    
for e in E:
    k1 = cmath.sqrt(2*m*e/(h**2))
    k2 = cmath.sqrt((2*m*(e-v0))/(h**2))
    k3 = k1

    try:
        M = matrix_method(k1,k2,k3,len(barriers))
        # print(M)

        r = -(M[1][0]/M[1][1])
        t = det(M)/M[1][1]

        R = r*np.conjugate(r)
        T = t*np.conjugate(t)

        Rs.append(R)
        Ts.append(T)

        total.append(T+R)
    except ZeroDivisionError:
         E = np.setdiff1d(E, np.array(e))

# plt.plot(Rs, E, label="Reflection")
plt.plot(Ts, E, label="Transmission")
# plt.plot(total, E)
plt.ylabel('Energy')
plt.title('Energy vs Transmission and Reflective Cofficients')
plt.xlabel('R, T')
plt.legend()
plt.show()