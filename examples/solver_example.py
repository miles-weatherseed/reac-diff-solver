import reac_diff_solver as rd
import numpy as np

# Define the initial conditions to be used by the solver. This is a function that takes two 2d numpy arrays which
# contain the x & y coordinates at each grid point and returns a list of two 2d numpy arrays (of the same shape) which
# contain the initial conditions for u and v.
def initial_conditions(X,Y):
    # uniform solution
    u = np.ones_like(X)
    v = np.zeros_like(X)

    # add perturbation in corner
    u[:20,:20] = 0.5
    v[:20,:20] = 0.25

    # add small amount of random noise to break the symmetry
    u += 0.1*np.random.standard_normal(X.shape)
    v += 0.1*np.random.standard_normal(X.shape)

    return [u,v]

# Define the reaction function to be used by the solver. This represents the reaction terms in the equation. The reaction
# function should take two 2d numpy arrays, which contain the values of u and v, and a list of parameters. It then returns
# a list of the values of the reaction terms f and g in the reaction-diffusion equation.
def Gray_Scott_reaction_terms(u,v, parameters = [1.0, 1.0]):
    k = parameters[0]
    F = parameters[1]
    return [-u*v**2 + F*(1-u), u*v**2 - (k+F)*v]

# Initialise the Solver. This creates a 256 x 256 grid running from x=0 to x=2.5 and y=0 to y=2.5 and calculates the
# initial conditions at each point on this grid.
solver = rd.Solver([0.0, 2.5], [0.0,2.5], 256, initial_conditions)

# Set the reaction terms.
solver.set_reactionFunction(Gray_Scott_reaction_terms)

# Set the timestep to be used for the integration.
solver.set_timeStepLength(10)

# Create a numpy array of times at which to find the solution. This creates an array with 100 timepoints between t = 0
# and t = 30000
t = np.linspace(0,30000,100)

# Run the Solver. The solve function takes the array of times and a list of the parameters to be used. The first two
# parameters are the diffusion coefficients of the u component and the v component respectively. The rest of the
# parameters are passed to the reaction function. In this case they represent the parameters k and F in the Gray-Scott
# equations. The final parameter tells the Solver to print out its progress occasionally.
u, v = solver.solve(t,[2E-5, 1E-5, 0.063, 0.032], True)

# Some other fun parameters to try are k = 0.050, F = 0.040 and k = 0.050, F = 0.018 (reduce the total time duration to
# 5000, and the timestep to 1)

# Create an animation of the solution and save it to a file called 'GrayScottEquationsAnimation'
rd.animate(u,10,"GrayScottEquations")