import numpy as np
import scipy.sparse
from ConjugateGradients import conjugate_gradients

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
        self.timeTrace = []
        self.dim = len(self.spatialLowerBounds)
        self.funcdim = 2
        self.funcparams = (1,)
        self.initialConditions = None

    def set_reactionFunction(self, function):
        self.reactionFunction = function

    def set_initialConditions(self, ic_function):
        self.initialCondtions = ic_function

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


    def make_arrays(self):
        """Make spatial mesh including final point"""

        self.tTrace = np.linspace(self.timeBounds[0], self.timeBounds[1], self.timeStepNumber)
        self.xTrace = np.linspace(self.xBounds[0], self.xBounds[1], self.xStepNumber)
        self.yTrace = np.linspace(self.yBounds[0], self.yBounds[1], self.yStepNumber)

        self.uSolution = np.zeros(self.timeStepNumber + 1, self.yStepNumber + 1, self.xStepNumber + 1)
        self.vSolution = np.zeros(self.timeStepNumber +1 , self.yStepNumber + 1, self.xStepNumber + 1)

    def set_initialconditions(self):
        self.uSolution[0,:,:] = self.initialConditions[0]
        self.vSolution[0,:,:] = self.initialConditions[1]

    def make_1dFdMatrix(self, N):
        """ Make FD Matrix for Periodic BCs"""
        nx = self.xStepNumber
        e = np.ones(nx + 1)
        diagonals = [e, -2 * e, e]
        offsets = [-1, 0, 1]
        Lx = (self.xStepLength)**2 * scipy.sparse.spdiags(diagonals, offsets, nx+1, nx+1,
                                 format='csr').tolil()
        # make periodic
        Lx[0, -1] = 1
        Lx[-1, 0] = 1

        ny = self.yStepNumber
        Ly = (self.yStepLength)**2 * scipy.sparse.spdiags(diagonals, offsets, ny + 1, ny + 1,
                                 format='csr').tolil()
        # make periodic
        Ly[0, -1] = 1
        Ly[-1, 0] = 1

        self.matrix = scipy.sparse.kron(np.eye(ny+1), Lx) + scipy.sparse.kron(Ly, np.eye(nx+1))



    def solve(self):
        uvec = self.uSolution[0,:,:].reshape(-1)
        vvec = self.vSolution[0,:,:].reshape(-1)
        I = np.eye(len(uvec))
        matrix = I - self.timeStepLength * self.matrix
        x0 = np.zeros(len(uvec))

        for i in range(0, self.timeStepNumber):
            uvec = self.uSolution[i,:,:].reshape(-1)
            vvec = self.vSolution[i,:,:].reshape(-1)
            uvec1 = conjugate_gradients(matrix, uvec + self.timeStepLength * self.reactionFunction[0](uvec, vvec), x0)
            vvec1 = conjugate_gradients(matrix, vvec + self.timeStepLength * self.reactionFunction[0](uvec, vvec), x0)

            self.uSolution[i+1,:,:] = uvec1.reshape(self.xStepNumber + 1, self.yStepNumber + 1)
            self.vSolution[i + 1, :, :] = vvec1.reshape(self.xStepNumber + 1, self.yStepNumber + 1)



