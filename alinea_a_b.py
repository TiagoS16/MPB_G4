import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from scipy.optimize import basinhopping


##BATCH##
def bl21_B(t, y, params):
    X, S, A, P = y
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, V0 = params
    u1 = 0.25 * (S / (0.3 + S))
    u2 = 0.55 * (S / (0.3 + S))
    u3 = 0.25 * (A / (0.4 + A))
    reac = [u1*X + u2*X + u3*X, -k1*u1*X - k2*u2*X, k3*u2*X - k4*u3*X, k11*u1*X] #X, S, A, P || usando as reacoes de crescimento
    return reac



# Initial conditions
X0= 7 #g/L
S0= 12 #g/L
A0= 0 #g/L
P0= 0 #g/L

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
V0= 9 #como é batch o volume não se altera


#lista com os valores iniciais fornecida a func
y0= [X0, S0, A0, P0]

#lista com os parametros fornecida a func
params= [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, V0]

t0= 0 #tempo inicial
t= 7 #tempo final
dt= 0.001 #intervalo de tempo entre reads


# Call the ODE solver
r = ode(bl21_B).set_integrator('lsoda', method='bdf', lband= 0) #lband é o limite inferior -- para nao haver valores negativos
r.set_initial_value(y0, t0).set_f_params(params)


#storing variables
T, x, s, a, p= [], [], [], [], []

while r.successful() and r.t < t:
    time= r.t + dt
    T.append(r.t)
    x.append(r.integrate(time)[0])
    s.append(r.integrate(time)[1])
    a.append(r.integrate(time)[2])
    p.append(r.integrate(time)[3])
    #print(time, r.integrate(time))


# using the storing variables to plot
T = np.array([0] + T)
x = np.array([y0[0]] + x)
s = np.array([y0[1]] + s)
a = np.array([y0[2]] + a)
p = np.array([y0[3]] + p)

'''
print('#'*40)
print(T)
print('#'*40)
print(x)
print('#'*40)
print(s)
print('#'*40)
print(a)
print('#'*40)
print(p)
print('#'*40)
'''

plt.plot(T, x, label='Biomassa', color='blue')
plt.plot(T, s, label='Substrato', color='red')
plt.plot(T, a, label='Acetato', color='green')
plt.plot(T, p, label='Proteína', color='pink')
plt.plot(T, [V0] * len(T), label='Volume (L)', color='Purple')
plt.legend(loc='best')
plt.xlabel('Tempo (h)')
plt.ylabel('Concentração (g/L)')
plt.xlim(0., 0.5)
plt.grid()
plt.show()