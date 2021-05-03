import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from scipy.optimize import basinhopping


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
    k4, u2, Ks3 = params

    u1 = 0.25 * (S / (0.3 + S))
    #u2 = 0.55 * (S / (0.3 + S))
    u3 = 0.25 * (A / (Ks3 + A))

    D = 0.7 / V
    #reac = [u1 * X + u2 * X + u3 * X - D * X, -k1 * u1 * X - k2 * u2 * X - D * S + D * Se, k3 * u2 * X - k4 * u3 * X - D * A, k11 * u1 * X - D * P, Fe]  # X, S, A, P, V || usando as reacoes de crescimento
    reac = [u1 * X + u2 * X + u3 * X - D * X, -4.412 * u1 * X - 22.22 * u2 * X - D * S + D * 350, 8.61 * u2 * X - k4 * u3 * X - D * A, 13.21 * u1 * X - D * P, 0.7]  # X, S, A, P, V || usando as reacoes de crescimento
    return reac


# Initial conditions
X0= 1 #g/L
S0= 0 #g/L
A0= 0 #g/L
P0= 0 #g/L
V= 3 #L

'''
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
'''

#lista com os valores iniciais fornecida a func
y0= [X0, S0, A0, P0, V]

#lista com os parametros fornecida a func
#params_old= [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, V0, Fe, Se]
k4= 9.846
u2 = 0.55 * (0 / (0.3 + 0))
Ks3= 0.4
params= [k4, u2, Ks3]

# Final time and step
t0= 0 #tempo inicial
t= 20 #tempo final
dt= 0.5 #intervalo de tempo entre reads


def estimate(params):

    """

    This will be our estimate function that works out as the calculation of the difference between the experimental
    and predicted values and can be used as the objective function

    :param params: parameters; array-like data structure (list, tuple, numpy array) for the ode
    :return: the error between measured and predicted data, i.e. difS + difX + difA + difV
    """

    # Consider using global to access and change global variables outside of this function
    # Otherwise, model, time, initial conditions and experimental data can be hardcoded here within the function
    # Nevertheless, use always the global Y (array to store the results) to re-write this variable with new results
    global model
    global t
    global t0
    #global tf
    global dados_exp
    global y0
    global Y
    global diferenca


    # ode
    # consider scipy.integrate.ode method, with the integrator lsoda and method bdf
    # you should apply the initial values and parameters
    # Call the ODE solver
    r = ode(model).set_integrator('lsoda', method='bdf', lband=0)  # lband é o limite inferior -- para nao haver valores negativos
    r.set_initial_value(y0, t0).set_f_params(params)

    # Using the global storing variable Y

    # Solve ode for the time step
    # see scipy docs on how to resolve the ode for the time steps
    # append the results to the Y storing variable

    # Y= np.array([1, 0, 0, 0, 3]) #storing with initial
    Y = [[1, 0, 0, 0, 3]]

    # T, x, s, a, p, v = [], [], [], [], [], []  # storing variables
    while r.successful() and r.t < t:
        time = r.t + dt
        Y.append(r.integrate(time).tolist())
        # T.append(r.t)
        # x.append(r.integrate(time)[0])
        # s.append(r.integrate(time)[1])
        # a.append(r.integrate(time)[2])
        # p.append(r.integrate(time)[3])
        # v.append(r.integrate(time)[4])
        # print(time, r.integrate(time))

    for i in range(len(Y)):
        Y[i].pop(3)

    for i in range(len(dados_exp)):
        dados_exp[i].pop(0)

    #Consider the metrics to calculate the error between experimental and predicted data
    diferenca= np.subtract(dados_exp,Y)

    #por aqui a formula para calcular o erro e por essa formula como output?
    return diferenca


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

# model = jm109
model = jm109

#t = timespan
#Final time and step
t1 = 20
dt = 0.5
#t= timespan(t0, t1, dt)


# dados_exp = pd.read_excel or pd.read_csv
dados_exp= pd.read_excel('dados_exp.xlsx').to_numpy().tolist()

# y0 = initial conditions

# Y = initialize the storing array. consider np.zeros()

# Y[0,:] = append the initial conditions to the results


######################################################################################
#Bounds
LB = [0, 0, 0]
UB = [4, 4, 4]
bounds = Bounds(LB, UB)

# Simulated Annealing
# consider scipy.optimize.basinhopping method, with the method BFGS, niter 200, seed 1 and the respective bounds.
# To perform some testing consider lowering the number of iterations to 10, as SA can take a while

# initial guess, that is the initial values for the parameters to be estimated. It can be those available in the pdf
x0 = [9.846, 0.55 * (0 / (0.3 + 0)), 0.4]

minimizer_kwargs = {"method": "BFGS"}

aserio= basinhopping(estimate, x0, minimizer_kwargs= minimizer_kwargs, niter= 200, accept_test= bounds, seed= 1)
#tentar= basinhopping(estimate, x0, minimizer_kwargs= minimizer_kwargs, niter= 10, accept_test= bounds)
#print(tentar)


print(estimate(params))