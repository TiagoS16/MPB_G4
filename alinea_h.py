import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sympy
from scipy.optimize import fsolve
from sympy import *


#sympy symbols
X, S, A, P = symbols('X S A P')

#Reactions
u1 = 0.25 * (S / (0.3 + S))
u2 = 0.55 * (S / (0.3 + S))
u3 = 0.25 * (A / (0.4 + A))

# Initial conditions
X0= 4 #g/L
S0= 0 #g/L
A0= 0 #g/L
P0= 0 #g/L
V= 8 #L

# Parameters
k1= 4.412
k2= 22.22
k3= 8.61
k4= 9.846
k5= 3.253
k6= 12.29
k7= 4.085
k8= 3.345
k9= 21.04
k10= 7.65
k11= 13.21
V0= 8 #o volume inicial não se altera
Se= 350 #concentração do substrato de entrada g/L
Fin= 0.7 #L/h || caudal de entrada || 350 g/L glucose
Fout= Fin #caudal de saida


difX= u1 * X + u2 * X + u3 * X - Fin * X
difS= - k1 * u1 * X - k2 * u2 * X + Fin * S0 - Fin * S
difA= k3 * u2 * X - k4 * u3 * X - Fin * A
difP= k11 * u1 * X - Fin * P
difV= Fin - Fout
todas= [difX, difS, difA, difP]


def derX(t=False):
    dX = [sympy.diff(difX, X), sympy.diff(difX, S), sympy.diff(difX, A), sympy.diff(difX, P)]
    if t == True:
        print('Derivada X X:')
        print(sympy.diff(difX, X))
        print()
        print('Derivada X S:')
        print(sympy.diff(difX, S))
        print()
        print('Derivada X A:')
        print(sympy.diff(difX, A))
        print()
        print('Derivada X P:')
        print(sympy.diff(difX, P))
    return dX


def derS(t=False):
    dS = [sympy.diff(difS, X), sympy.diff(difS, S), sympy.diff(difS, A), sympy.diff(difS, P)]
    if t == True:
        print('Derivada S X:')
        print(sympy.diff(difS, X))
        print()
        print('Derivada S S:')
        print(sympy.diff(difS, S))
        print()
        print('Derivada S A:')
        print(sympy.diff(difS, A))
        print()
        print('Derivada S P:')
        print(sympy.diff(difS, P))
    return dS


def derA(t=False):
    dA = [sympy.diff(difA, X), sympy.diff(difA, S), sympy.diff(difA, A), sympy.diff(difA, P)]
    if t == True:
        print('Derivada A X:')
        print(sympy.diff(difA, X))
        print()
        print('Derivada A S:')
        print(sympy.diff(difA, S))
        print()
        print('Derivada A A:')
        print(sympy.diff(difA, A))
        print()
        print('Derivada A P:')
        print(sympy.diff(difA, P))
    return dA


def derP(t=False):
    dP = [sympy.diff(difP, X), sympy.diff(difP, S), sympy.diff(difP, A), sympy.diff(difP, P)]
    if t == True:
        print('Derivada P X:')
        print(sympy.diff(difP, X))
        print()
        print('Derivada P S:')
        print(sympy.diff(difP, S))
        print()
        print('Derivada P A:')
        print(sympy.diff(difP, A))
        print()
        print('Derivada P P:')
        print(sympy.diff(difP, P))
    return dP

derX()
print('='*66)
derS()
print('='*66)
derA()
print('='*66)
derP()


def fun(X):
    S= 0
    A= 0
    P= 0
    # Reactions
    u1 = 0.25 * (S / (0.3 + S))
    u2 = 0.55 * (S / (0.3 + S))
    u3 = 0.25 * (A / (0.4 + A))

    # Initial conditions
    S0 = 0  # g/L

    # Parameters
    k1 = 4.412
    k2 = 22.22
    k3 = 8.61
    k4 = 9.846
    k11 = 13.21
    Fin = 0.7

    return [(u1 * X[0] + u2 * X[0] + u3 * X[0] - Fin * X[0]),
            (- k1 * u1 * X[1] - k2 * u2 * X[1] + Fin * S0 - Fin * S),
            (k3 * u2 * X[2] - k4 * u3 * X[2] - Fin * A),
            (k11 * u1 * X[3] - Fin * P)]


root = fsolve(fun, [1, 1, 1, 1])
print(root)
print(np.isclose(fun(root), [0.0, 0.0, 0.0, 0.0]))
