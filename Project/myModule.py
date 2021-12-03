"""
myModule.py
We create this module to have all classes and functions we created during solving problem of our memoir
"""

# RK4 function definition
from copy import deepcopy

import numpy as np
from matplotlib import pyplot as plt


def RK4(f, I, y0, n, sol_only: bool = False):
    """
    This function allows to solve Ordinary differential Equation using RK4 method
    :param f: a single or a system of ordinary differential equation as a function
    :param I: interval which solution will be searched
    :param y0: initial values as a vector
    :param n: number of interval subdivision
    :param sol_only: True (default: False) if only the solution and not with derivative solution will be returned
    :return: a vector of subdivided interval and the solution of the ODE
    """
    h = (I[1] - I[0]) / n
    t = np.linspace(I[0], I[1], n + 1)
    m = y0.size
    y = np.zeros((n + 1, m))
    y[0, :] = y0
    for i in range(n):
        k1 = h * f(t[i], y[i, :])
        k2 = h * f(t[i] + h / 2, y[i, :] + k1 / 2)
        k3 = h * f(t[i] + h / 2, y[i, :] + k2 / 2)
        k4 = h * f(t[i] + h, y[i, :] + k3)

        y[i + 1, :] = y[i, :] + (k1 + 2 * (k2 + k3) + k4) / 6
        if sol_only:
            return t, y[:, 0]
    return t, y


def Euler(f, I, y0, n,):
    h = (I[1] - I[0]) / n
    t = np.linspace(I[0], I[1], n + 1)
    m = y0.size
    y = np.zeros((n + 1, m))
    y[0, :] = y0


def pop_dyn_daudi(t, pop, rates):
    rates = deepcopy(rates)
    pop = deepcopy(pop)
    t = deepcopy(t)

    x = pop[0]  # Maize population
    y = pop[1]  # larvae population
    z = pop[2]  # noctudae population
    w = pop[3]  # eggs population

    an = rates[0]  # $\alpha$ or $\eta$
    l = rates[1]  # $\lambda$
    ga = rates[2]  # $\gamma$
    d = rates[3]  # $\delta$
    rh = rates[4]  # $\rho$
    my = rates[5]  # $\mu_y$
    mz = rates[6]  # $\mu_z$
    mw = rates[7]  # $\mu_w$
    e = rates[8]  # $e$

    f = lambda x: -l * x #* (x / 500 - 1)
    g = lambda x: an * x

    dx = f(x) - g(x) * y
    dy = e * g(x) * y + ga * w - (d + my) * y
    dz = d * y - mz * z
    dw = rh * z - (ga + mw) * w

    return np.array([dx, dy, dz, dw])


def pop_dyn_myMod(t, pop, rates):
    rates = deepcopy(rates)
    pop = deepcopy(pop)
    t = deepcopy(t)

    x = pop[0]  # Maize population
    y = pop[1]  # larvae population
    z = pop[2]  # adult population
    w = pop[3]  # eggs population

    an = rates[0]  # destruction of plant per larvae over time
    l = rates[1]  # maize mortality rate due to climatic condition
    ga = rates[2]  # egg -> larvae
    d = rates[3]  # larvae -> adult
    rh = rates[4]  # fertility rate
    my = rates[5]  # larvae mortality rate
    mz = rates[6]  # adult mortality rate
    mw = rates[7]  # egg mortality rate
    e = rates[8]  # larvae growth or survival rate per plant consumption
    tau = rates[9]  # resistance rate of the plant over time
    sig = rates[10]  # migration rate of adult
    s = rates[11]  # threshold of maize decreasing

    myf = lambda x: -l * x * (x / s - 1)  # growth function
    myg = lambda x: an * x  # / (tau*t + 1)    # functional response

    dx = (myf(x) - myg(x) * y)  # * np.exp(-tau * t) #/ (tau * t + 1)  # maize population dynamic over time
    dy = -(e * myg(x) + d) * y + ga * w - my * y  # larvae population dynamic
    dz = (e * myg(x) + d) * y - mz * z  # adult population dynamic
    dw = rh * z - (ga + mw) * w  # egg population dynamic

    return np.array([dx, dy, dz, dw])


class solver:
    def __init__(self, rate, s0, I, sub=None):
        self.rate = deepcopy(rate)
        self.s0 = deepcopy(s0)
        self.I = deepcopy(I)
        self.sub = deepcopy(I[-1]) * 100 if sub is None else deepcopy(sub)
        self.t = None
        self.sol = None
        self.title = ["Maize", "Caterpillar", "Adult moth", "Egg"]

    def myModSol(self):
        myf = lambda t, x: pop_dyn_myMod(t, x, self.rate)
        self.modSol(myf)

    def daudiModSol(self):
        daudf = lambda t, x: pop_dyn_daudi(t, x, self.rate)
        self.modSol(daudf)

    def modSol(self, f):
        self.t, self.sol = RK4(f, self.I, self.s0, self.sub)

    def set_input(self, rate, s0, I, sub=None):
        self.__init__(rate, s0, I, sub)

    def printSol(self, fig):
        for i, ax in enumerate(fig.axes):
            ax.plot(self.t, self.sol[:, i], label="z0 = {}".format(self.s0[2]))
            ax.set_title(self.title[i])
        fig.tight_layout()
        return fig
