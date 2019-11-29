import numpy as np
import scipy.sparse
from src.conjgrad import conjugate_gradients


class Solver:
        # Constructor
    def __init__(self, xBounds, yBounds, gridSize, initial_condition_function):
        self.initial_condition_function = initial_condition_function
        self.set_grid(xBounds,yBounds,gridSize)
        self.initialConditions_u = self.initial_condition_function(self.X, self.Y)[0]
        self.initialConditions_v = self.initial_condition_function(self.X, self.Y)[1]

        self.reactionFunction = lambda u, v, parameters: [0,0]
        self.timeStepLength = 10

    def set_grid(self, xBounds, yBounds, gridSize):
        self.xBounds = xBounds
        self.yBounds = yBounds
        self.gridSize = gridSize
        self.xStepLength = (xBounds[1] - xBounds[0])/(gridSize) # +1 improves solution but inconsistent ???
        self.yStepLength = (yBounds[1] - yBounds[0])/(gridSize)
        self.x = np.linspace(self.xBounds[0], self.xBounds[1], gridSize+1)[:-1]
        self.y = np.linspace(self.yBounds[0], self.yBounds[1], gridSize+1)[:-1]
        self.X, self.Y = np.meshgrid(self.x, self.y) # create mesh
        self.initialConditions = self.initial_condition_function(self.X, self.Y) # calculate initial conditions on new mesh

    def set_reactionFunction(self, function):
        self.reactionFunction = function

    def set_initialConditions(self, initial_condition_function):
        self.initial_condition_function = initial_condition_function
        self.initialConditions_u = self.initial_condition_function(self.X, self.Y)[0]
        self.initialConditions_v = self.initial_condition_function(self.X, self.Y)[1]

    def set_timeStepLength(self, length):
        self.timeStepLength = length

    def _create_arrays(self, times):
        """Create arrays to store solution at given times"""

        self.uSolution = np.zeros((len(times), self.gridSize, self.gridSize))
        self.vSolution = np.zeros((len(times), self.gridSize, self.gridSize))


    def _create_fdmatrix(self):
        """ Make FD Matrix for Periodic BCs"""
        n = self.gridSize
        e = np.ones(n)
        diagonals = [e, -2 * e, e]
        offsets = [-1, 0, 1]
        L = scipy.sparse.spdiags(diagonals, offsets, n, n, format='csr').tolil()

        # periodic boundary conditions
        L[0, -1] = 1
        L[-1, 0] = 1

        I = scipy.sparse.eye(n)

        self.laplacian = scipy.sparse.kron(I, self.xStepLength**-2 * L) + scipy.sparse.kron(self.yStepLength**-2 * L, I)


    def solve(self, times, parameters):
        # first two parameters are diffusion coefficients
        self.diff_u = parameters[0]
        self.diff_v = parameters[1]

        self._create_arrays(times)
        self._create_fdmatrix()

        self.uSolution[0] = self.initialConditions_u
        self.vSolution[0] = self.initialConditions_v

        uvec = self.uSolution[0].reshape(-1)
        vvec = self.vSolution[0].reshape(-1)

        I = scipy.sparse.eye(self.gridSize**2)
        matrix_u = (I - self.timeStepLength * self.diff_u * self.laplacian)
        matrix_v = (I - self.timeStepLength * self.diff_v * self.laplacian)

        t = times[0]
        for i in range(1,len(times)):
            j = 1
            print(times[i])
            while t < times[i]:
                uvec = conjugate_gradients(matrix_u, uvec + self.timeStepLength * self.reactionFunction(uvec, vvec, parameters[2:])[0], uvec)[0]
                vvec = conjugate_gradients(matrix_v, vvec + self.timeStepLength * self.reactionFunction(uvec, vvec, parameters[2:])[1], vvec)[0]
                t = times[i-1] + j*self.timeStepLength
                j += 1

            self.uSolution[i] = uvec.reshape(self.gridSize, self.gridSize)
            self.vSolution[i] = vvec.reshape(self.gridSize, self.gridSize)

        return [self.uSolution, self.vSolution]
