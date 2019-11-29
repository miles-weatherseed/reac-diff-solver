from src.solver import Solver
import numpy as np
from src.animation import animate

def initial_conditions(X,Y):
    u = np.ones_like(X)
    v = np.zeros_like(X)
    u[:20,:20] = 0.5
    v[:20,:20] = 0.25
    u += 0.1*np.random.standard_normal(X.shape)
    v += 0.1*np.random.standard_normal(X.shape)

    return [u,v]


def reaction_function(u,v, parameters = [1.0, 1.0]):
    F = parameters[1]
    k = parameters[0]
    return [-u*v**2 + F*(1-u), u*v**2 - (k+F)*v]


solver = Solver([0.0, 2.5], [0.0,2.5], 256, initial_conditions)

solver.set_reactionFunction(reaction_function)
t = np.linspace(0,30000,100)
u, v = solver.solve(t,[2E-5, 1E-5, 0.063, 0.032])
T,Y,X = np.meshgrid(t,solver.y,solver.x, indexing = 'ij')

animate(u,10)