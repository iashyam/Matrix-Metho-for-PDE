import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.linalg import det

v0 = 5
m = 1
h = 1
a = 5
E = np.linspace(0,10, 200)

i = 1j #defining iota it is not good for me to work with j

def scatterinMatrix(k1, k2):
    d = np.zeros([2,2], dtype=complex)
    if k2 == 0:
        k2 = e-10
    d[0][0] = (1+k1/k2)
    d[0][1] = (1-k1/k2)
    d[1][1] = (1+k1/k2)
    d[1][0] = (1-k1/k2)

    return (1/2)*d

def progogationMatrix(a, k):

    return np.array([[cmath.exp(i*k*a),0],[0 , cmath.exp(-i*k*a)]], dtype=complex)

def matrixMethod(E):
    k1 = cmath.sqrt(2*m*E/(h**2))
    k2 = cmath.sqrt((2*m*(E-v0))/(h**2))
    k3 = k1

    d21 = scatterinMatrix(k1,k2)
    d32 = scatterinMatrix(k2,k3)

    p = progogationMatrix(a, k2)

    M =  np.matmul(d32, np.matmul(p,d21))

    r = -(M[1][0]/M[1][1])
    t = det(M)/M[1][1]

    R = r*np.conjugate(r)
    T = t*np.conjugate(t)

    return R, T

r = []
t = []
total = []

for e in E:
    r.append(matrixMethod(e)[0])
    t.append(matrixMethod(e)[1])

for y in range(len(r)):
    total.append(r[y] + t[y])



plt.plot(E, t)
# plt.plot(E,total)
plt.show()