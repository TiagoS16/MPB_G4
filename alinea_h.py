import pandas as pd
import numpy as np
import numdifftools as nd
import matplotlib.pyplot as plt
import sympy as sym
from scipy.optimize import fsolve
from sympy import *


#sympy symbols
X, S, A, P = symbols('X,S,A,P')#,umax1,umax2,umax3,Ks1,Ks2,Ks3,Fin,Fout,S0,k1,k2,k3,k4,k11

#Reactions
# u1 = umax1 * (S / (Ks1 + S))
# u2 = umax2 * (S / (Ks2 + S))
# u3 = umax3 * (A / (Ks3 + A))
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
    dX = [sym.diff(difX, X), sym.diff(difX, S), sym.diff(difX, A), sym.diff(difX, P)]
    if t == True:
        print('Derivada X X:')
        print(sym.diff(difX, X))
        print()
        print('Derivada X S:')
        print(sym.diff(difX, S))
        print()
        print('Derivada X A:')
        print(sym.diff(difX, A))
        print()
        print('Derivada X P:')
        print(sym.diff(difX, P))
    return dX


def derS(t=False):
    dS = [sym.diff(difS, X), sym.diff(difS, S), sym.diff(difS, A), sym.diff(difS, P)]
    if t == True:
        print('Derivada S X:')
        print(sym.diff(difS, X))
        print()
        print('Derivada S S:')
        print(sym.diff(difS, S))
        print()
        print('Derivada S A:')
        print(sym.diff(difS, A))
        print()
        print('Derivada S P:')
        print(sym.diff(difS, P))
    return dS


def derA(t=False):
    dA = [sym.diff(difA, X), sym.diff(difA, S), sym.diff(difA, A), sym.diff(difA, P)]
    if t == True:
        print('Derivada A X:')
        print(sym.diff(difA, X))
        print()
        print('Derivada A S:')
        print(sym.diff(difA, S))
        print()
        print('Derivada A A:')
        print(sym.diff(difA, A))
        print()
        print('Derivada A P:')
        print(sym.diff(difA, P))
    return dA


def derP(t=False):
    dP = [sym.diff(difP, X), sym.diff(difP, S), sym.diff(difP, A), sym.diff(difP, P)]
    if t == True:
        print('Derivada P X:')
        print(sym.diff(difP, X))
        print()
        print('Derivada P S:')
        print(sym.diff(difP, S))
        print()
        print('Derivada P A:')
        print(sym.diff(difP, A))
        print()
        print('Derivada P P:')
        print(sym.diff(difP, P))
    return dP

derX()
#print('='*66)
derS()
#print('='*66)
derA()
#print('='*66)
derP()


def fun(x):
    X, S, A, P = x
    # Reactions
    u1 = 0.25 * (S / (0.3 + S))
    u2 = 0.55 * (S / (0.3 + S))
    u3 = 0.25 * (A / (0.4 + A))

    # Initial conditions
    S0 = 12 # g/L

    # Parameters
    k1 = 4.412
    k2 = 22.22
    k3 = 8.61
    k4 = 9.846
    k11 = 13.21
    Fin = 0.7

    return [(u1 * X + u2 * X + u3 * X - Fin * X),
            (- k1 * u1 * X - k2 * u2 * X + Fin * S0 - Fin * S),
            (k3 * u2 * X - k4 * u3 * X - Fin * A),
            (k11 * u1 * X - Fin * P)]


root = fsolve(fun, [0, 12, 0, 0])
print(root)
print(np.isclose(fun(root), [0.0, 0.0, 0.0, 0.0]))
print()
root1 = fsolve(fun, [7, 12, 0, 0])
print(root1)
print(np.isclose(fun(root1), [0.0, 0.0, 0.0, 0.0]))
print()


# difX= u1 * X + u2 * X + u3 * X - Fin * X
# difS= - k1 * u1 * X - k2 * u2 * X + Fin * S0 - Fin * S
# difA= k3 * u2 * X - k4 * u3 * X - Fin * A
# difP= k11 * u1 * X - Fin * P
# difV= Fin - Fout


def Jacobian(v_str, f_list):
    global J
    vars = sym.symbols(v_str)
    f = sym.sympify(f_list)
    J = sym.zeros(len(f),len(vars))
    for i, fi in enumerate(f):
        for j, s in enumerate(vars):
            J[i,j] = sym.diff(fi, s)
    J = np.array(J).tolist()
    return J
Jacobian('X G A P', ['u1*X + u2*X + u3*X - Fin*X',' - k1*u1*X - k2*u2*X + Fin*S0 - G*Fin', 'k3*u2* X - k4*u3*X - Fin*A', 'k11*u1*X - Fin*P'])

#ponto 1
#X1, G1, A1, P1 = 0, 12, 0, 0

#u1 = 0.2439
#u2 = 0.53659
#u3 = 0

JJ = [[0.08049, 0, 0, 0], [-12.9991, -0.7, 0, 0], [4.62, 0, -0.7, 0], [3.22192, 0, 0, -0.7]]
JJJ = np.array(JJ)
jjjj = np.linalg.det(JJJ)
print('O determinante do ponto 1 é {}'.format(jjjj))
jjjjj = np.trace(JJJ)
print('O traço do ponto 1 é {}'.format(jjjjj))


#ponto 2
#X2, S2, A2, P2 = -1.37289662*(10**-9), 1.20000000*(10**1), -1.06074905*(10**-8), -6.31708043*(10**-9)
#u1=0.2439
#u2=0.5638
#u3=0.1524

LL = [[0.2611, 0, 0, 0], [-13.6037, -0.7, 0, 0], [3.3538, 0, -0.7, 0], [3.22192, 0, 0, -0.7]]
LLL = np.array(LL)
llll = np.linalg.det(LLL)
print('O determinante do ponto 2 é {}'.format(llll))
lllll = np.trace(LLL)
print('O traço do ponto 2 é {}'.format(lllll))


#são os dois pontos sela