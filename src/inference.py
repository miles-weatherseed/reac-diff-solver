from scipy import optimize
from solver import Solver

class Inference:
    def __init__(self, data, times):
        '''
        Takes a 2 x n numpy array of generated values of u and v
        :param data: the generated data
        :param times: the corresponding times
        '''
        self.data = data
        self.times = times

    def error_func(self, output):
        '''
        returns L2 error between values generated using current set of parameters and the original data
        :param output:
        :return:
        '''
        tot_error = 0
        for i in range(len(self.times)):
            tot_error += np.dot(self.data[:, i] - output[:,i],  self.data[:, i] - output[:,i])
        return tot_error

    def cost_func(self, parametervalues):
        solver = Solver()
        self.output = solver.solve(self.times, parametervalues)
        return self.error_func(self.output)

    def fit_params(self, x0):
        res = optimize.minimize(self.cost_func, x0, method='Nelder-Mead')
        return res.x