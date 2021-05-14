import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import sympy


def bl21_FB(t, y, params):
    X, S, A, P, V = y
    k1, k2, k3, k5, k6, k7, k8, k9, k10, k11, V0, Fe, Se = params
    u1 = 0.25 * (S / (0.3 + S))
    u2 = 0.55 * (S / (0.3 + S))
    u3 = 0.25 * (A / (0.4 + A))
    D= Fe / V
    reac = [u1*X + u2*X + u3*X - D*X, -k1*u1*X - k2*u2*X - D*S + D*Se, k3*u2*X - 9.846*u3*X - D*A, k11*u1*X - D*P, Fe] #X, S, A, P, V || usando as reacoes de crescimento
    return reac


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
k5= 3.253
k6= 12.29
k7= 4.085
k8= 3.345
k9= 21.04
k10= 7.65
k11= 13.21
V0= 8 #o volume inicial não se altera
Fe= 0.7 #L/h || caudal de entrada || 350 g/L glucose
Se= 350 #concentração do substrato de entrada g/L

#lista com os valores iniciais fornecida a func
y0= [X0, S0, A0, P0, V]

#lista com os parametros fornecida a func
params= [k1, k2, k3, k5, k6, k7, k8, k9, k10, k11, V0, Fe, Se]

t0= 0 #tempo inicial
t= 20 #tempo final
dt= 0.001 #intervalo de tempo entre reads

#umax1 = 0.25
#umax3 = 0.25
#Ks1 = 0.3
#Ks2 = 0.3

# Consider using sympy.symbols to create algebric variables to be used on the derivatives (X, S, k1, ks1, ...)
#(X,S,A,P,V,k4,umax2,Ks3) = sympy.symbols('X,S,A,P,V,k4,umax2,ks3', real=True)
(X,S,A,P,V,k1,k2,k4,umax1,umax2,umax3,Ks1,Ks2,Ks3,Fe,Se) = sympy.symbols('X,S,A,P,V,k1,k2,k4,umax1,umax2,umax3,Ks1,Ks2,Ks3,Fe,Se', real=True)


u1 = umax1 * (S / (Ks1 + S))
u2 = umax2 * (S / (Ks2 + S))
u3 = umax3 * (A / (Ks3 + A))

#ordinary differential equations DX and DS
dxdt = u1*X + u2*X + u3*X - Fe/V*X
dsdt = -k1*u1*X - k2*u2*X - Fe/V*S + Fe/V*Se

#derivadas parciais
dX_k4 = sympy.diff(dxdt, k4)
dX_umax2 = sympy.diff(dxdt, umax2)
dX_Ks3 = sympy.diff(dxdt, Ks3)
dS_k4 = sympy.diff(dsdt, k4)
dS_umax2 = sympy.diff(dsdt, umax2)
dS_Ks3 = sympy.diff(dsdt, Ks3)



#funções lambda -- lambdify
dX_k4 = sympy.lambdify((), dX_k4, 'numpy')  # derivada parcial = 0
dX_umax2 = sympy.lambdify((X, S, Ks2), dX_umax2, 'numpy')
dX_Ks3 = sympy.lambdify((X, A, umax3, Ks3), dX_Ks3, 'numpy')
dS_k4 = sympy.lambdify((), dS_k4, 'numpy')
dS_umax2 = sympy.lambdify((X, S, k2, Ks2), dS_umax2, 'numpy')
dS_Ks3 = sympy.lambdify((), dS_Ks3, 'numpy')


# Call the ODE solver
r = ode(bl21_FB).set_integrator('lsoda', method='bdf', lband=0) #lband é o limite inferior -- para nao haver valores negativos
r.set_initial_value(y0, t0).set_f_params(params)

#storing variables
T, x, s, a, p, v= [], [], [], [], [], []

while r.successful() and r.t < t:
    time= r.t + dt
    T.append(r.t)
    x.append(r.integrate(time)[0])
    s.append(r.integrate(time)[1])
    a.append(r.integrate(time)[2])
    p.append(r.integrate(time)[3])
    v.append(r.integrate(time)[4])
    #print(time, r.integrate(time))


T = np.array([0] + T)
x = np.array([y0[0]] + x)
s = np.array([y0[1]] + s)
a = np.array([y0[2]] + a)
p = np.array([y0[3]] + p)
v = np.array([y0[4]] + v)


X = x
S = s
A = a
k2 = 22.22
k4 = 9.846  # restaurar os valores originais 
umax2 = 0.55  # restaurar os valores originais
umax3 = 0.25
Ks3 = 0.4  # restaurar os valores originais
Ks2 = 0.3


# create the dx and ds arrays executing the lambda function for the X, S, k4, umax2, ks3 values
dx_k4 = np.zeros(20001)
dx_umax2 = dX_umax2(X, S, Ks2)
dx_Ks3 = dX_Ks3(X, A, umax3, Ks3)
ds_k4 = np.zeros(20001)
ds_umax2 = dS_umax2(X, S, k2, Ks2)
ds_Ks3 = np.zeros(20001)


### plot X
fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, dx_k4, linewidth=2, linestyle='solid', label='X_K4', color='blue')
ax.set_xlabel('Tempo (h)')
ax.set_title('Sensibilidade de X a K4 - BL21 (Fed-Batch)')
ax.legend()
plt.grid()
plt.show()

fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, dx_umax2, '--', linewidth=2, label='X_umax2', color='red')
ax.set_xlabel('Tempo (h)')
ax.set_title('Sensibilidade de X a umax2 - BL21 (Fed-Batch)')
ax.legend()
plt.grid()
plt.show()

fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, dx_Ks3, '--', linewidth=2, label='X_Ks3', color='orange')
ax.set_xlabel('Tempo (h)')
ax.set_title('Sensibilidade de X a Ks3 - BL21 (Fed-Batch)')
ax.legend()
plt.grid()
plt.show()

fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, dx_k4, linewidth=2, linestyle='solid', label='X_K4', color='blue')
ax.plot(T, dx_umax2, '--', linewidth=2, label='X_umax2', color='red')
ax.plot(T, dx_Ks3, '--', linewidth=2, label='X_Ks3', color='orange')
ax.set_xlabel('Tempo (h)')
ax.set_title('Sensibilidades de X - BL21 (Fed-Batch)')
ax.legend()
plt.grid()
plt.show()

### plot S
fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, ds_k4, linewidth=2, linestyle='solid', label='S_K4', color='blue')
ax.set_xlabel('Tempo (h)')
ax.set_title('Sensibilidade de S a K4 - BL21 (Fed-Batch)')
ax.legend()
plt.grid()
plt.show()

fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, ds_umax2, '--', linewidth=2, label='S_umax2', color='red')
ax.set_xlabel('Tempo (h)')
ax.set_title('Sensibilidade de S a umax2 - BL21 (Fed-Batch)')
ax.legend()
plt.grid()
plt.show()

fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, ds_Ks3, '--', linewidth=2, label='S_Ks3', color='orange')
ax.set_xlabel('Tempo (h)')
ax.set_title('Sensibilidade de S a Ks3 - BL21 (Fed-Batch)')
ax.legend()
plt.grid()
plt.show()

fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, ds_k4, linewidth=2, linestyle='solid', label='S_K4', color='blue')
ax.plot(T, ds_umax2, '--', linewidth=2, label='S_umax2', color='red')
ax.plot(T, ds_Ks3, '--', linewidth=2, label='S_Ks3', color='orange')
ax.set_xlabel('Tempo (h)')
ax.set_title('Sensibilidades de S - BL21 (Fed-Batch)')
ax.legend()
plt.grid()
plt.show()