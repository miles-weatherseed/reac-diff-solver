from reac_diff_solver.solver import Solver
import numpy as np
from reac_diff_solver.animation import animate

# Example code to obtain solution for reaction diffusion equation of Fizthugh-Nagumo type
# epsilon * u_t = du * (u_xx + u_yy) -u^3 + 3 * u - v
# v_t = dv * (v_xx + v_yy) + u - delta * v
# where du, dv are diffusion coefficients

# Define Range of x and y values
upperbound = 100
xbounds = [0, upperbound]
ybounds = [0,upperbound]

# Define Initial Condition Function to obtain spiralling asymptotic behaviour
# Quadrants taking values 0, 1, 0, -1
def initial_conditions_spiral(X,Y):
    u = 0.5 * (X >= upperbound/2) + 0.5 * (Y >= upperbound/2) - 0.5 * (Y < upperbound/2) - 0.5 * (X < upperbound/2)
    v = 0.5 * (X >= upperbound/2) + 0.5 * (Y >= upperbound/2) - 0.5 * (Y < upperbound/2) - 0.5 * (X < upperbound/2)
    return [u, v]

# Define Random Initial Condition Function to obtain wave like asymptotic behaviour
def initial_conditions_random(X, Y):
    u = 0.1*np.random.standard_normal(X.shape)
    v = 0.1*np.random.standard_normal(X.shape)
    return [u, v]

# Define parameters for Fitzhugh-Nagumo reaction term
eps = 0.1       #epsilon
delta = 0.001   #delta
du = 1.0        # (diffusion coefficient)
dv = 1.0        # (diffusion coefficient)

# Rewrite diffusion coefficients to match the form of solver
Du = 1/eps
Dv = dv

# Define Reaction Function
def fitzhugh_naguomo(u, v, parameters):
    e = parameters[0]
    d = parameters[1]
    return [(3*u - u**3 -v)/e, u - d * v]


# Initialise Solver
# Set x, y bounds and grid size
# Set initial condition function to obtain spiralling behavior
solver = Solver(xbounds, ybounds, 256, initial_conditions_spiral)

# Set reaction terms
solver.set_reactionFunction(fitzhugh_naguomo)

# Set time step lengths for backward Euler integration scheme
solver.set_timeStepLength(0.01)

# Define array of times at which to record solution
t = np.linspace(0,30,100)

# Run solve method of solver.
# Input time array of times to record solution and list of parameters for the system.
# First two parameters are diffusion coefficients, remaining parameters are passed to diffusion function
# Return solution arrays for u, v
u, v = solver.solve(t, [Du, Dv, eps, delta])
T,Y,X = np.meshgrid(t,solver.y,solver.x, indexing = 'ij')

# Produce animation of function and save to file as 'FitzhughNagumoAnimation'
animate(u,10, "FitzhughNagumoAnimation")