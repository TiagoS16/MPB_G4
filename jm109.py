import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from scipy.optimize import basinhopping


def jm109(t, Y, params):
    '''
    X, S, A, P, V = y
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, V0, Fe, Se = params
    u1 = 0.25 * (S / (0.3 + S))
    u2 = 0.55 * (S / (0.3 + S))
    u3 = 0.25 * (A / (0.4 + A))
    D = Fe / V
    reac = [u1 * X + u2 * X + u3 * X - D * X, -k1 * u1 * X - k2 * u2 * X - D * S + D * Se,
            k3 * u2 * X - k4 * u3 * X - D * A, k11 * u1 * X - D * P,
            Fe]  # X, S, A, P, V || usando as reacoes de crescimento
    return reac
    '''
    pass


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
    global dados_exp
    global y0
    global Y

    # timespan (initial, final, step)
    # Initial time
    t0 = 0
    t = 20
    dt = 0.5

    # Initial conditions
    X0 =   # g/L
    S0 =   # g/L
    A0 =   # g/L
    P0 =   # g/L
    V =   # L

    # ode
    # consider scipy.integrate.ode method, with the integrator lsoda and method bdf
    # you should apply the initial values and parameters

    # Final time and step
    # t1 =
    # dt =

    # Using the global storing variable Y

    # Solve ode for the time step
    # see scipy docs on how to resolve the ode for the time steps
    # append the results to the Y storing variable

    # Consider the metrics to calculate the error between experimental and predicted data

    #return
    pass


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

# t = timespan
t0= 0
t= 20
dt= 0.5


# dados_exp = pd.read_excel or pd.read_csv
dados_exp= pd.read_excel('dados_exp.xlsx').to_numpy()

# y0 = initial conditions
y0= [X0, S0, A0, P0, V]

# Y = initialize the storing array. consider np.zeros()

# Y[0,:] = append the initial conditions to the results

#Bounds
LB = [0, 0, 0]
UB = [4, 4, 4]
bounds = Bounds(LB, UB)

# Simulated Annealing
# consider scipy.optimize.basinhopping method, with the method BFGS, niter 200, seed 1 and the respective bounds.
# To perform some testing consider lowering the number of iterations to 10, as SA can take a while

# initial guess, that is the initial values for the parameters to be estimated. It can be those available in the pdf
x0 = []

# Consider using matplotlib for plotting.



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



estimar= basinhopping()