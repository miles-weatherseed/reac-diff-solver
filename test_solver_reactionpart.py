from solver import Solver
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def initial_conditions(X,Y):
    return [np.ones_like(X), np.zeros_like(X)]


def reaction_function(u,v, parameters = [1.0]):
    # parameters[0]: angular velocity
    w = parameters[0]
    return [w*v, -w*u]


solver = Solver([0.0,1.0], [0.0,1.0], 20, initial_conditions)
solver.set_reactionFunction(reaction_function)
t = np.linspace(0,2*np.pi,100)
u, v = solver.solve(t,[.5,1.0])

T,Y,X = np.meshgrid(t,solver.y,solver.x, indexing = 'ij')

plt.plot(t, u[:,10,10], label = 'u')
plt.plot(t, v[:,10,10], label = 'v')
plt.legend()
plt.show()