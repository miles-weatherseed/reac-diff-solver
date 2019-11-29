from src.solver import Solver
import numpy as np
import matplotlib.pyplot as plt


def test_1d_diffusion():
    mode = 1.0
    D1 = 0.05
    D2 = 0.0
    t = np.linspace(0, 1.0, 11)

    def exact(x, y, t):
        return np.sin(4 * np.pi * mode * x) * np.exp(-D1 * t * (4 * np.pi * mode) ** 2)

    def initial_conditions(x, y):
        return [exact(x, y, 0), 1.0]

    solver = Solver([0.0, 1.0], [0.0, 1.0], 200, initial_conditions)
    solver.set_timeStepLength(0.001)

    u, v = solver.solve(t, [D1, D2, 0])

    T, Y, X = np.meshgrid(t, solver.y, solver.x, indexing='ij')
    solution = exact(X, Y, T)

    plt.figure()
    for i in range(0, len(t)):
        plt.plot(solver.x, solution[i, 20, :])
    plt.show()

    plt.figure()
    for i in range(0, len(t)):
        plt.plot(solver.x, u[i, 20, :])

    plt.show()


    error = solution - u
    print(np.max(error))

def test_1d_reaction():
    mode = 1.0
    D1 = 0.05
    D2 = 0.0
    t = np.linspace(0, 1.0, 11)
    k1 = 1.0
    k2 = 2.0

    def exact(x, y, t):
        return [1.0 * np.exp(-k1*t), 2.0*np.exp(-k2*t)]

    def initial_conditions(x, y):
        return exact(x, y, 0)

    solver = Solver([0.0, 1.0], [0.0, 1.0], 200, initial_conditions)
    solver.set_timeStepLength(0.001)

    def reaction_function(u, v, K):
        return [-K[0] * u, - K[1] * v]

    solver.set_reactionFunction(reaction_function)

    u, v = solver.solve(t, [D1, D2, k1, k2])

    T, Y, X = np.meshgrid(t, solver.y, solver.x, indexing='ij')
    solution = exact(X, Y, T)

    plt.figure()
    plt.plot(t, u[:, 20, 20])
    plt.plot(t, v[:, 20, 20])

    plt.show()



#def reaction_function(u,v):
 #   return [0, v*(2-v)]





"""
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
"""

if __name__ == '__main__':
    test_1d_reaction()
    #test_1d_diffusion()


