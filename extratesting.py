from solver import Solver
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

mode = 1.0
D1 = 0.05
D2 = 0.0
def zero_mode_solution(x, y, t):
    return np.sin(4*np.pi*mode*x)*np.exp(-D1*t*(4*np.pi*mode)**2)

def initial_conditions(x, y):
    return [zero_mode_solution(x, y, 0), 1.0]

def reaction_function(u,v):
    return [0, v*(2-v)]

solver = Solver([0.0, 1.0], [0.0, 1.0], 200, initial_conditions)
t = np.linspace(0,0.1,11)
solver.set_reactionFunction(reaction_function)
u,v = solver.solve(t,[D1,D2])

T,Y,X = np.meshgrid(t,solver.y,solver.x, indexing = 'ij')
solution = zero_mode_solution(X, Y, T)

plt.figure()
plt.plot(solver.x, solution[1,0,:])
plt.plot(solver.x, u[1,100,:])
error = solution[1,0,:] - u[1,100,:]
print(np.max(error))
plt.show()

print(solver.diff_v)

plt.figure()
plt.plot(solver.x, len(solver.x)*[np.exp(-0.1)])
plt.plot(solver.x, v[-1,100,:])
error = solution[1,0,:] - v[1,100,:]
print(np.min(v[1,100,:]))
plt.show()

plt.figure()
plt.plot(t, np.exp(-t))
plt.plot(t, v[:,100,100])
#error = solution[1,0,:] - v[1,100,:]
#print(np.min(v[1,100,:]))
plt.show()

m = 1.0
n = 2.0
def exact_2d(x, y, t):
    return np.sin(2 * m * np.pi * x) * np.sin(2 * n * np.pi * y) * np.exp(-2 * (n**2 + m**2) * np.pi ** 2 * t)

def initial_conditions_2d(x,y):
    return [exact_2d(x, y, 0)] * 2

solver = Solver([0.0,1.0], [0.0,1.0], 200, initial_conditions_2d)
t = np.array([0.0,0.1])
u,v = solver.solve(t,[1.0,0.5])




T,Y,X = np.meshgrid(t,solver.y,solver.x, indexing = 'ij')
solution = exact_2d(X, Y, T)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(solver.X,solver.Y,v[1])
ax.plot_surface(solver.X,solver.Y,solution[1])
plt.show()
