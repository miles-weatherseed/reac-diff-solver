from scipy import optimize
from solver import Solver
import numpy as np

'''
We can take a time series which we believe to have been obtained from a reaction diffusion PDE system and run this
inference class on it. The output will be the parameters that minimize some cost function defined as L2 error.
'''

class Inference(Solver):
    def __init__(self, u_data, v_data, times):
        '''
        Initializes our inference class by receiving the observed grids of u and v values and the times in the time
        series. We set the x and y boundaries and the initial conditions separately.

        :param u_data: observed values of u at each point in the grid at each time in the time series
        :type u_data: n x gridsize x gridsize numpy array
        :param u_data: observed values of v at each point in the grid at each time in the time series
        :type u_data: n x gridsize x gridsize numpy array
        :param times: array of n evenly-spaced times
        :type times: 1-dimensional numpy array
        '''

        self.u_data = u_data
        self.v_data = v_data
        self.times = times

        self.gridsize = u_data.shape[1] #infer the gridsize from this dimension

    def set_model(self, xBounds, yBounds, initial_conditions_function):
        '''
        Sets up the reaction diffusion model we wish to consider. We determine the x and y boundaries for consideration
        and provide some initial conditions for our system.

        :param xBounds: x-range of the problem
        :type u_data: list of 2 floats
        :param yBounds: y-range of the problem
        :type u_data: list of 2 floats
        :param initial_condition_function: calculates the values of u and v at t=0 at each gridpoint
        :type initial_condition_function: function that takes two 2d numpy arrays and returns a list of two 2d numpy arrays
        '''

        self.solver = Solver(xBounds, yBounds, self.gridsize, initial_conditions_function)

    def set_reaction_function(self, function):
        '''
        This is optional. We can choose to provide a nonlinear term to our system.
        :param function: calculates the value of the reaction terms at the given u and v
        :type function: function that takes two numpy arrays (containing values of u and v) and a list of parameters
        and returns a list of two numpy arrays
        '''

        self.solver.set_reactionFunction(function)

    def error_func(self):
        '''
        Returns L2 error between output using current proposed set of parameters and the observed data
        :return tot_error: the total sum of L2 errors at each timestep between the true u,v values and the solved values
        '''

        tot_error = 0
        for i in range(len(self.times)):
            u_diff_mat = self.u_data[i, :, :] - self.u_output[i, :, :]
            v_diff_mat = self.v_data[i, :, :] - self.v_output[i, :, :]
            tot_error += np.sum(np.square(u_diff_mat)) + np.sum(np.square(v_diff_mat))
        return tot_error

    def cost_func(self, parametervalues):
        '''
        Takes newly proposed set of parameter values and returns the L2 error between the values of u,v given by the
        solver with these values and the values of u,v from the input
        :param parametervalues: the proposed set of parameter values at this point in the optimization
        :return: the L2 error with these parameters
        '''

        self.u_output, self.v_output = self.solver.solve(self.times, parametervalues)
        return self.error_func()

    def fit_params(self, x0):
        '''
        The master function to call. Fits parameters from some initial estimate by minimizing the cost function using
        the Nelder-Mead method.
        :param x0: the initial estimate of parameters
        :return best_fit_parameters: a vector of the fitted parameters
        '''
        
        res = optimize.minimize(self.cost_func, x0, method='Nelder-Mead')
        best_fit_parameters = res.x
        return best_fit_parameters