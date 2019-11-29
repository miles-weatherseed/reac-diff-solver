from solver import Solver
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def returns_zero(u,v):
    return 0

def initial_conditions(X,Y):
    return [np.sin(4*np.pi*X) * np.sin(2*np.pi*Y), np.sin(4*np.pi*X) * np.sin(2*np.pi*Y)]

solver = Solver([0.0,1.0], [0.0,1.0], 200, initial_conditions)
t = np.array([0.0,0.1])
u,v = solver.solve(t,[.5,1.0])

T,Y,X = np.meshgrid(t,solver.y,solver.x, indexing = 'ij')
solution = np.sin(4*np.pi*X) * np.sin(2*np.pi*Y) * np.exp(-20*np.pi**2 * T)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(solver.X,solver.Y,v[1])
ax.plot_surface(solver.X,solver.Y,solution[1])
plt.show()


print(np.max(v[1] - solution[1]))