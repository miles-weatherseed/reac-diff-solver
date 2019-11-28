import numpy as np
import scipy.sparse
from ConjugateGradients import conjugate_gradients
import matplotlib.pyplot as plt


class Solver:

        # Constructor
    def __init__(self):
        self.timeBounds = [0.0, 1.0]
        self.xBounds = [0.0, 1.0]
        self.yBounds = [0.0, 1.0]
        self.reactionFunction = None
        self.timeStepLength = 0.001
        self.timeStepNumber = None
        self.xStepLength = 0.01
        self.yStepLength = 0.01
        self.xStepNumber = None
        self.yStepNumber = None
        self.solutionTrace = []
        self.tTrace = []
        self.funcdim = 2
        self.funcparams = (1,)
        self.initialConditions = None

    def set_reactionFunction(self, function):
        self.reactionFunction = function

    def set_initialConditions(self, initial_conditions):
        self.initialConditions = initial_conditions

    def set_timeBounds(self, bounds):
        self.timeBounds = bounds

    def set_xBounds(self, bounds):
        self.xBounds = bounds

    def set_yBounds(self, bounds):
        self.yBounds = bounds

    def set_timeStepLength(self, length):
        self.timeStepLength = length

    def set_xStepLength(self, length):
        self.xStepLength = length

    def set_yStepLength(self, length):
        self.yStepLength = length

    def get_step_numbers(self):
        """
        Determine number of steps for finite differences and updates corresponding step sizes
        """

        self.timeStepNumber = int(np.ceil((self.timeBounds[1] - self.timeBounds[0]) / self.timeStepLength))
        self.timeStepLength = (self.timeBounds[1] - self.timeBounds[0]) / self.timeStepNumber

        self.xStepNumber = int(np.ceil((self.xBounds[1] - self.xBounds[0]) / self.xStepLength))
        self.xStepLength = (self.xBounds[1] - self.xBounds[0]) / self.xStepNumber

        self.yStepNumber = int(np.ceil((self.yBounds[1] - self.yBounds[0]) / self.yStepLength))
        self.yStepLength = (self.yBounds[1] - self.yBounds[0]) / self.yStepNumber


    def create_arrays(self):
        """Make spatial mesh including final point"""


        self.tTrace = np.linspace(self.timeBounds[0], self.timeBounds[1], self.timeStepNumber)
        self.xTrace = np.linspace(self.xBounds[0], self.xBounds[1], self.xStepNumber)
        self.yTrace = np.linspace(self.yBounds[0], self.yBounds[1], self.yStepNumber)

        self.uSolution = np.zeros((self.timeStepNumber, self.yStepNumber, self.xStepNumber))
        self.vSolution = np.zeros((self.timeStepNumber, self.yStepNumber, self.xStepNumber))


    def create_fdmatrix(self):
        """ Make FD Matrix for Periodic BCs"""
        nx = self.xStepNumber
        e = np.ones(nx)
        diagonals = [e, -2 * e, e]
        offsets = [-1, 0, 1]
        Lx = (self.xStepLength)**-2 * scipy.sparse.spdiags(diagonals, offsets, nx, nx,
                                 format='csr').tolil()
        # make periodic
        Lx[0, -1] = (self.xStepLength)**-2
        Lx[-1, 0] = (self.xStepLength)**-2

        ny = self.yStepNumber
        Ly = (self.yStepLength)**-2 * scipy.sparse.spdiags(diagonals, offsets, ny, ny,
                                 format='csr').tolil()
        # make periodic
        Ly[0, -1] = (self.yStepLength)**-2
        Ly[-1, 0] = (self.yStepLength)**-2

        self.matrix = scipy.sparse.kron(scipy.sparse.eye(ny), Lx) + scipy.sparse.kron(Ly, scipy.sparse.eye(nx))


    def solve(self):
        self.get_step_numbers()
        self.create_arrays()
        self.create_fdmatrix()

        self.uSolution[0, :, :] = self.initialConditions[0]
        self.vSolution[0, :, :] = self.initialConditions[1]

        uvec = self.uSolution[0,:,:].reshape(-1)
        I = scipy.sparse.eye(len(uvec))
        matrix = I - self.timeStepLength * self.matrix
        x0 = np.zeros(len(uvec))

        for i in range(1, self.timeStepNumber):
            uvec = self.uSolution[i-1,:,:].reshape(-1)
            vvec = self.vSolution[i-1,:,:].reshape(-1)
            uvec1 = conjugate_gradients(matrix, uvec + self.timeStepLength * self.reactionFunction[0](uvec, vvec), x0)[0]
            vvec1 = conjugate_gradients(matrix, vvec + self.timeStepLength * self.reactionFunction[1](uvec, vvec), x0)[0]

            self.uSolution[i ,: , :] = uvec1.reshape(self.xStepNumber, self.yStepNumber)
            self.vSolution[i ,: , :] = vvec1.reshape(self.xStepNumber, self.yStepNumber)

