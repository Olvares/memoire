# ================= IMPORTATION =================
from copy import deepcopy
# %matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import myModule as myM

#%%
leg = ['Maize', 'Caterpillar', 'Adult', 'Eggs']

# ================= DAUDI's MODEL SOLUTION =================
# initialisation
k = 500
y0 = 0
z0 = np.array([15, 30, 45, 60])
w0 = 0

r = 0.0417
d = 0.071
ga = 0.071
my = 0.0071
mz = 0.115
mw = 0.04
l = 0.015
e = 1.6
an = 0.000154

rate = np.array([an, l, ga, d, r, my, mz, mw, e])

t0 = 0
t1 = 63
T = 160
I1 = [t0, T]

i = 0
s0 = np.array([k, y0, z0[i], w0])
sol = myM.solver(rate, s0, I1)
# t, x = sol.daudiModSol()
# print(p)

plt.figure(1)
for i in range(z0.size):
    s0 = np.array([k, y0, z0[i], w0])
    sol.set_input(rate, s0, I1)
    sol.daudiModSol()
    sol.printSol(plt)
    # for j in range(4):
    #     plt.subplot(2, 2, j+1)
    #     plt.plot(t, x[:, j])
    # plt.title(leg[j])
    # plt.legend(["z0 = {}".format(z0[i])])

plt.show()

#%%
# =================+++ MY MODEL SOLUTION ++=================
# initialisation
k = 300
y0 = 0
z0 = np.array([15, 30, 45, 60])
w0 = 0

r = 0.10#0.0417     # fertility rate
d = 0.071           # larvae -> adult
ga = 0.071          # egg -> larvae
my = 0.0071         # larvae mortality rate
mz = 0.115          # adult mortality rate
mw = 0.04           # egg mortality rate
l = 0#.0015 #0.015   # maize mortality rate due to climatic condition
e = 1.6             # larvae growth or survival rate per plant consumption
an = 0.000154       # destruction of plant per larvae over time
tau = 0.1  # resistance rate of the plant over time
sig = 0.01  # migration rate of adult
s = k / 3  # threshold

rate = np.array([an, l, ga, d, r, my, mz, mw, e, tau, sig, s])

# t0 = 0
# t1 = 63
# T = 160
#
# I1 = [t0,t1]

plt.figure(2)
for i in range(z0.size):
    s0 = np.array([k, y0, z0[i], w0])
    sol.set_input(rate, s0, I1)
    sol.myModSol()
    sol.printSol(plt)
    # for j in range(4):
    #     plt.subplot(2, 2, j+1)
    #     plt.plot(t, x[:, j])
    # plt.title(leg[j])
    # plt.legend(["z0 = {}".format(z0[i])])
plt.show()
plt.close("all")