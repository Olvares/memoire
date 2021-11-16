"""
myModule.py
We create this module to have all classes and functions we created during solving problem of our memoir
"""

# RK4 function definition
import numpy as np


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


def pop_dyn_daudi(t, pop, rates=None):
    if rates is None:
        rates = []
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

    f = lambda x: l * x * (1 - x / 400)
    g = lambda x: an * x
    dx = f(x) - g(x) * y  # -(an*y + l) * x
    dy = e * g(x) * y + ga * w - (d + my) * y
    dz = d * y - mz * z
    dw = rh * z - (ga + mw) * w

    return np.array([dx, dy, dz, dw])


def pop_dyn_myMod(t, pop, rates=None):
    if rates is None:
        rates = []
    x = pop[0]  # Maize population
    y = pop[1]  # larvae population
    z = pop[2]  # adult population
    w = pop[3]  # eggs population

    an = rates[0]       # destruction of plant per larvae over time
    l = rates[1]        # maize mortality rate due to climatic condition
    ga = rates[2]       # egg -> larvae
    d = rates[3]        # larvae -> adult
    rh = rates[4]       # fertility rate
    my = rates[5]       # larvae mortality rate
    mz = rates[6]       # adult mortality rate
    mw = rates[7]       # egg mortality rate
    e = rates[8]        # larvae growth or survival rate per plant consumption
    tau = rates[9]      # resistance rate of the plant over time
    sig = rates[10]     # migration rate of adult
    s = rates[11]       # threshold of maize decreasing

    f = lambda x: -l * x * (x / s - 1)    # growth function
    g = lambda x: an * x #/ (tau*t + 1)    # functional response

    dx = f(x) - g(x) * y #/ (tau*t + 1)                       # maize population dynamic over time
    dy = -(e * g(x) + d) * y + ga * w #- my * y   # larvae population dynamic
    dz = (e * g(x) + d) * y - mz * z                         # adult population dynamic
    dw = rh * z - (ga + mw) * w                 # egg population dynamic

    return np.array([dx, dy, dz, dw])
