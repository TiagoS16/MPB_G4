import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from scipy.optimize import basinhopping

dados_exp= pd.read_excel('dados_exp.xlsx').to_numpy().tolist()
for i in range(len(dados_exp)): #retira a coluna do tempo(T) para ter uma matriz igual à Y dos estimados
    dados_exp[i].pop(0)


def bl21_FB(t, y, params):
    X, S, A, P, V = y
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, V0, Fe, Se, umax2, Ks3 = params
    u1 = 0.25 * (S / (0.3 + S))
    u2 = umax2 * (S / (0.3 + S))
    u3 = 0.25 * (A / (Ks3 + A))
    D= Fe / V
    reac = [u1*X + u2*X + u3*X - D*X, -k1*u1*X - k2*u2*X - D*S + D*Se, k3*u2*X - k4*u3*X - D*A, k11*u1*X - D*P, Fe] #X, S, A, P, V || usando as reacoes de crescimento
    return reac


def run_ode(modelo, param):
    global y0
    r = ode(modelo).set_integrator('lsoda', method='bdf', lband=0, nsteps=5000)  # lband é o limite inferior -- para nao haver valores negativos
    r.set_initial_value(y0, t0).set_f_params(param)

    Y = [[1, 0, 0, 0, 3]]  # variavel Y com os dados iniciais

    while r.successful() and r.t < t:
        time = r.t + dt
        Y.append(r.integrate(time).tolist())

    for i in range(len(Y)):  # retira a coluna da proteina(P) para ter uma matriz igual à dos dados experimentais
        Y[i].pop(3)
    return Y

modelo = bl21_FB

# Initial conditions
X0 = 1 #g/L  BL21=4  dados_exp=1
S0 = 0 #g/L
A0 = 0 #g/L
P0 = 0 #g/L
V = 3 #L  BL21=8  dados_exp=3

#lista com os valores iniciais fornecida a func
y0 = [X0, S0, A0, P0, V]

# Final time and step
t0= 0 #tempo inicial
t= 20 #tempo final
dt= 0.5 #intervalo de tempo entre reads

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
V0= 3 #o volume inicial não se altera
Fe= 0.7 #L/h || caudal de entrada || 350 g/L glucose
Se= 350 #concentração do substrato de entrada g/L
umax2 = 0.55
Ks3 = 0.4

#lista com os parametros fornecida a func
param= [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, V0, Fe, Se, umax2, Ks3]

inicial = run_ode(modelo, param)

Yx, Ys, Ya, Yv = [], [], [], []
DEx, DEs, DEa, DEv = [], [], [], []
T = [0]

for i in range(40): #criar a lista com os tempos para fazer os graficos
    T.append(T[i]+0.5)

for i in range(len(inicial)): #separar as colunas da matriz dos estimados para obter os valores para fazer o grafico
    Yx.append(inicial[i][0])
    Ys.append(inicial[i][1])
    Ya.append(inicial[i][2])
    Yv.append(inicial[i][3])

for i in range(len(dados_exp)): #separar as colunas da matriz dos dados experimentais para obter os valores para fazer o grafico
    DEx.append(dados_exp[i][0])
    DEs.append(dados_exp[i][1])
    DEa.append(dados_exp[i][2])
    DEv.append(dados_exp[i][3])


#grafico conjunto inicial
fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, Yx, linewidth=2, label='Biomassa Modelo BL21', color='blue')
ax.plot(T, DEx, 'o-', markersize=4, linewidth=1, label='Biomassa Dados_exp', color='blue')
ax.plot(T, Ys, linewidth=2, label='Substrato Modelo BL21', color='red')
ax.plot(T, DEs, 'o-', markersize=4, linewidth=1, label='Substrato Dados_exp', color='red')
ax.plot(T, Ya, linewidth=2, label='Acetato Modelo BL21', color='green')
ax.plot(T, DEa, 'o-', markersize=4, linewidth=1, label='Acetato Dados_exp', color='green')
ax.plot(T, Yv, linewidth=2, label='Volume Modelo BL21 (L)', color='purple')
ax.plot(T, DEv, 'o-', markersize=4, linewidth=1, label='Volume Dados_exp (L)', color='purple')
ax.set_xlabel('Tempo (h)')
ax.set_ylabel('Concentração (g/L)')
ax.set_title('Modelo BL21 vs Dados_exp JM109 (Fed-Batch)')
ax.legend(loc='best')
plt.grid()
plt.show()


def jm109(t, y, params):
    '''
    This will be the model for the strain JM109 which is similar to the BL21, but it should have slight modifications
    :param t: time; This argument should not be altered
    :param Y: initial conditions; array-like data structure (list, tuple, numpy array)
    :param params: parameters; array-like data structure (list, tuple, numpy array) - NOTE THAT THESE ARGUMENT MUST
    CONTAIN ONLY AND ONLY THOSE PARAMETERS TO BE ESTIMATED. The remain parameters should be hardcoded within the
    function
    :return: K * phi - (D * variables) + zeros; note that numpy.dot() is the standard for matrices multiplication
    '''

    X, S, A, P, V = y
    k4, umax2, Ks3 = params

    u1 = 0.25 * (S / (0.3 + S))
    u2 = umax2 * (S / (0.3 + S))
    u3 = 0.25 * (A / (Ks3 + A))

    D = 0.7 / V
    reac = [u1 * X + u2 * X + u3 * X - D * X, -4.412 * u1 * X - 22.22 * u2 * X - D * S + D * 350, 8.61 * u2 * X - k4 * u3 * X - D * A, 13.21 * u1 * X - D * P, 0.7]  # X, S, A, P, V || usando as reacoes de crescimento
    return reac


def estimate(params):
    """
    This will be our estimate function that works out as the calculation of the difference between the experimental
    and predicted values and can be used as the objective function

    :param params: parameters; array-like data structure (list, tuple, numpy array) for the ode
    :return: the error between measured and predicted data, i.e. difS + difX + difA + difV
    """

    global model
    global t
    global t0
    global dt
    global dados_exp
    global y0
    global Y
    global soma1
    global Y_legit

    r = ode(model).set_integrator('lsoda', method='bdf', lband=0, nsteps=5000)  # lband é o limite inferior -- para nao haver valores negativos
    r.set_initial_value(y0, t0).set_f_params(params)

    Y = [[1, 0, 0, 0, 3]] #variavel Y com os dados iniciais

    while r.successful() and r.t < t:
        time = r.t + dt
        Y.append(r.integrate(time).tolist())

    for i in range(len(Y)): #retira a coluna da proteina(P) para ter uma matriz igual à dos dados experimentais
        Y[i].pop(3)

    if len(Y) == len(dados_exp): #a função só vai executar isto quando os tamanhos das 2 matrizes forem iguais
        #isto é necessário porque a ODE estava a terminar com uma matriz com menos linhas e não permitia a execução do calculo do erro
        #o anterior acontecia porque a ODE não estava a conseguir integrar com sucesso (while r.successful no chunk acima)
        #para evitar isto fizemos com que o basinhopping use os valores anteriores da ODE caso a atual execução da mesma dê o tal erro
        Y_legit= Y #guardar numa variavel diferente para no final conseguir usar a matriz dos estimados sem as falhas faladas anteriormente
        #diferenca= np.subtract(Y, dados_exp)
        diferenca = np.subtract(Y_legit, dados_exp) #diferença entre os valores
        po = np.power(diferenca, 2) #diferença ao quadrado
        soma= po.sum(axis=0) #soma dos quadrados das diferenças
        soma1= soma.sum() #soma das somas dos quadrados das diferenças
    return soma1 #soma dos quadrados das diferenças


# Bounds
# Consider using the following class for setting the Simulated Annealing bounds
class Bounds(object):

    def __init__(self, LB=None, UB=None):

        if LB is None:
            LB = [0, 0, 0]

        if UB is None:
            UB = [4, 4, 4]

        self.lower_bound = np.array(LB)
        self.upper_bound = np.array(UB)

    def __call__(self, **kwargs):

        x = kwargs["x_new"]

        tmax = bool(np.all(x <= self.upper_bound))
        tmin = bool(np.all(x >= self.lower_bound))

        return tmax and tmin


#Bounds
LB = [0, 0, 0]
UB = [4, 4, 4]
bounds = Bounds(LB, UB)

model = jm109

# initial guess, that is the initial values for the parameters to be estimated. It can be those available in the pdf
x0 = [9.846, 0.55, 0.4]

minimizer_kwargs = {"method": "BFGS"} #method BFGS

for_real = basinhopping(estimate, x0, minimizer_kwargs=minimizer_kwargs, niter=200, accept_test=bounds, seed=1)
#tentar = basinhopping(estimate, x0, minimizer_kwargs=minimizer_kwargs, accept_test=bounds, niter=1, seed=1) #niter_success para a otimização caso o mínimo se mantenha igual em n iterações sucessivas
param_est = for_real.x
print(for_real)
print('Os mínimos encontrados são {}.'.format(param_est))

#redefinir os parametros para os minimos estimados
k4 = param_est[0]
umax2 = param_est[1]
Ks3 = param_est[2]

#lista com os parametros fornecida a func
param = [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, V0, Fe, Se, umax2, Ks3]

modelo = bl21_FB

final = run_ode(modelo, param)

Yx, Ys, Ya, Yv = [], [], [], []

for i in range(len(final)): #separar as colunas da matriz dos estimados para obter os valores para fazer o grafico
    Yx.append(final[i][0])
    Ys.append(final[i][1])
    Ya.append(final[i][2])
    Yv.append(final[i][3])


#grafico conjunto final
fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, Yx, linewidth=2, label='Biomassa Estimado JM109', color='blue')
ax.plot(T, DEx, 'o-', markersize=4, linewidth=1, label='Biomassa Dados_exp', color='blue')
ax.plot(T, Ys, linewidth=2, label='Substrato Estimado JM109', color='red')
ax.plot(T, DEs, 'o-', markersize=4, linewidth=1, label='Substrato Dados_exp', color='red')
ax.plot(T, Ya, linewidth=2, label='Acetato Estimado JM109', color='green')
ax.plot(T, DEa, 'o-', markersize=4, linewidth=1, label='Acetato Dados_exp', color='green')
ax.plot(T, Yv, linewidth=2, label='Volume Estimado JM109 (L)', color='purple')
ax.plot(T, DEv, 'o-', markersize=4, linewidth=1, label='Volume Dados_exp (L)', color='purple')
ax.set_xlabel('Tempo (h)')
ax.set_ylabel('Concentração (g/L)')
ax.set_title('Modelo Estimado JM109 vs Dados_exp JM109 (Fed-Batch)')
ax.legend(loc='best')
plt.grid()
plt.show()


#OLD GRAPHS
#plot do grafico do estimado
# fig, ax = plt.subplots()
# fig.set_figheight(8)
# fig.set_figwidth(10)
# ax.plot(T, Yx, label='Biomassa', color='blue')
# ax.plot(T, Ys, label='Substrato', color='red')
# ax.plot(T, Ya, label='Acetato', color='green')
# ax.plot(T, Yv, label='Volume (L)', color='Purple')
# ax.legend(loc='best')
# ax.set_xlabel('Tempo (h)')
# ax.set_ylabel('Concentração (g/L)')
# ax.set_title('Modelo JM109 (Fed-Batch)')
# plt.grid()
# plt.show()
#
#
# #plot do grafico dos dados experimentais
# fig, ax = plt.subplots()
# fig.set_figheight(8)
# fig.set_figwidth(10)
# ax.plot(T, DEx, label='Biomassa', color='blue')
# ax.plot(T, DEs, label='Substrato', color='red')
# ax.plot(T, DEa, label='Acetato', color='green')
# ax.plot(T, DEv, label='Volume (L)', color='purple')
# ax.legend(loc='best')
# ax.set_xlabel('Tempo (h)')
# ax.set_ylabel('Concentração (g/L)')
# ax.set_title('Dados Experimentais JM109 (Fed-Batch)')
# plt.grid()
# plt.show()