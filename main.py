import cmath
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import det


# THR INITIAL CONDITIONS
v0 = 5 #height of the potential box is 5
m = 1 #taking mass to be 1
h = 1 #taking plank constant to be 1
a = 1 #widht of potential box is 1
E = np.linspace(0, 10, 200)

Rs = []  # list of reflecing and trasmattting constant
Ts = []
total = []  # list of r+t

i = 1j  # defining iota it is not good for me to work with j

#scattering matrix depands in k1 and k2
def scatterinMatrix(k1, k2):
    d = np.zeros([2, 2], dtype=complex)
    if k2 == 0: #if k32 is 0 we can take in so small
        k2 = e-10
    d[0][0] = (1+k1/k2)
    d[0][1] = (1-k1/k2)
    d[1][1] = (1+k1/k2)
    d[1][0] = (1-k1/k2)

    return (1/2)*d


def progogationMatrix(a, k):

    return np.array([[cmath.exp(i*k*a), 0], [0, cmath.exp(-i*k*a)]], dtype=complex)


for e in E:
    k1 = cmath.sqrt(2*m*e/(h**2))
    k2 = cmath.sqrt((2*m*(e-v0))/(h**2))
    k3 = k1

    d21 = scatterinMatrix(k1, k2)
    d32 = scatterinMatrix(k2, k3)

    p = progogationMatrix(a, k2)

    M = np.matmul(d32, np.matmul(p, d21))

    r = -(M[1][0]/M[1][1])
    t = det(M)/M[1][1]

    R = r*np.conjugate(r)
    T = t*np.conjugate(t)

    Rs.append(R)
    Ts.append(T)
    total.append(T+R)

plt.plot(E, Ts, label='Transmission')
plt.plot(E, Rs, label="Reflection")
plt.plot(E, total)
plt.xlabel('Energy')
plt.title('Energy vs Transmission and Reflective Cofficients')
plt.ylabel('R, T')
plt.legend()
plt.show()