import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode


def bl21_FB(t, y, params):
    X, S, A, P, V = y
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, V0, Fe, Se = params
    u1 = 0.25 * (S / (0.3 + S))
    u2 = 0.55 * (S / (0.3 + S))
    u3 = 0.25 * (A / (0.4 + A))
    D= Fe / V
    reac = [u1*X + u2*X + u3*X - D*X,
            -k1*u1*X - k2*u2*X - D*S + D*Se,
            k3*u2*X - k4*u3*X - D*A,
            k11*u1*X - D*P,
            Fe] #X, S, A, P, V || usando as reacoes de crescimento
    return reac

def bl21_FB_redux(t, y, params):
    # remover k4 porque sensibilidade nula
    X, S, A, P, V = y
    k1, k2, k3, k5, k6, k7, k8, k9, k10, k11, V0, Fe, Se = params
    u1 = 0.25 * (S / (0.3 + S))
    u2 = 0.55 * (S / (0.3 + S))
    u3 = 0.25 * (A / (0.4 + A))
    D= Fe / V
    reac = [u1*X + u2*X + u3*X - D*X,
            -k1*u1*X - k2*u2*X - D*S + D*Se,
            k3*u2*X - 1*u3*X - D*A,
            k11*u1*X - D*P,
            Fe] #X, S, A, P, V || usando as reacoes de crescimento
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
k4= 9.846
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
params= [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, V0, Fe, Se]

t0= 0 #tempo inicial
t= 20 #tempo final
dt= 0.001 #intervalo de tempo entre reads

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


params= [k1, k2, k3, k5, k6, k7, k8, k9, k10, k11, V0, Fe, Se]
#ODE solver for redux
r = ode(bl21_FB_redux).set_integrator('lsoda', method='bdf', lband=0) #lband é o limite inferior -- para nao haver valores negativos
r.set_initial_value(y0, t0).set_f_params(params)

#storing variables
T_redux, x_redux, s_redux, a_redux, p_redux, v_redux= [], [], [], [], [], []

while r.successful() and r.t < t:
    time= r.t + dt
    T_redux.append(r.t)
    x_redux.append(r.integrate(time)[0])
    s_redux.append(r.integrate(time)[1])
    a_redux.append(r.integrate(time)[2])
    p_redux.append(r.integrate(time)[3])
    v_redux.append(r.integrate(time)[4])
    #print(time, r.integrate(time))


T_redux = np.array([0] + T_redux)
x_redux = np.array([y0[0]] + x_redux)
s_redux = np.array([y0[1]] + s_redux)
a_redux = np.array([y0[2]] + a_redux)
p_redux = np.array([y0[3]] + p_redux)
v_redux = np.array([y0[4]] + v_redux)

#Plot
fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T_redux, x_redux, linewidth=2, label='Biomassa_BL21_mod_reduzido', color='blue')
ax.plot(T_redux, s_redux, linewidth=2, label='Substrato_BL21_mod_reduzido', color='red')
ax.plot(T_redux, a_redux, linewidth=2, label='Acetato_BL21_mod_reduzido', color='green')
ax.plot(T_redux, p_redux, linewidth=2, label='Proteína Recombinada_BL21_mod_reduzido', color='pink')
ax.plot(T_redux, v_redux, linewidth=2, label='Volume_BL21_mod_reduzido (L)', color='purple')
ax.plot(T, x, '--', linewidth=2, label='Biomassa_BL21', color='blue')
ax.plot(T, s, '--', linewidth=2, label='Substrato_BL21', color='red')
ax.plot(T, a, '--', linewidth=2, label='Acetato_BL21', color='green')
ax.plot(T, p, '--', linewidth=2, label='Proteína Recombinada_BL21', color='pink')
ax.plot(T, v, '--', linewidth=2, label='Volume_BL21 (L)', color='purple')
ax.set_xlabel('Tempo (h)')
ax.set_ylabel('Concentrações (g/L)')
ax.set_title('Model BL21 (Fed-Batch)  VS Model BL21 (Fed-Batch) Reduzido')
ax.legend()
plt.grid()
plt.show()