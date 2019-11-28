from solver import Solver
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def returns_zero(u,v):
    return 0

solver = Solver()
solver.set_timeStepLength(0.0001)
solver.get_step_numbers()
solver.create_arrays()

x = solver.xTrace
y = solver.yTrace
X,Y = np.meshgrid(x,y)

initial_conditions = np.sin(4*np.pi*X) * np.sin(2*np.pi*Y)

solver.set_initialConditions([initial_conditions,initial_conditions])
solver.set_reactionFunction([returns_zero,returns_zero])
solver.solve()

t = solver.tTrace
T,Y,X = np.meshgrid(t,y,x, indexing = 'ij')
solution = np.sin(4*np.pi*X) * np.sin(2*np.pi*Y) * np.exp(-40*np.pi**2 * T)


X,Y = np.meshgrid(x,y)


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X,Y,solver.uSolution[250])
ax.plot_surface(X,Y,solution[250])
plt.show()