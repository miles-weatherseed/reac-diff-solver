from scipy import optimize
from solver import Solver
import numpy as np

'''
We can take a time series which we believe to have been obtained from a reaction diffusion PDE system and run this
inference class on it. The output will be the parameters
'''

class Inference(Solver):
    def __init__(self, u_data, v_data, times):
        '''
        Takes an n dimensional vector of times with corresponding observed data for u and v in (n, x, x) arrays
        where x is the dimension of the mesh. Initalizes these in the solver.
        :param data: the observed data we wish to infer the parameters from
        :param times: the corresponding times
        '''
        self.u_data = u_data
        self.v_data = v_data
        self.times = times

    def set_model(self, xBounds, yBounds, gridsize, initial_conditions_function):
        self.solver = Solver(xBounds, yBounds, gridsize, initial_conditions_function)

    def set_reaction_function(self, function):
        self.solver.set_reactionFunction(function)

    def error_func(self):
        '''
        Returns L2 error between output using current proposed set of parameters and the observed data
        :param u_output:
        :param v_output:
        :return:
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