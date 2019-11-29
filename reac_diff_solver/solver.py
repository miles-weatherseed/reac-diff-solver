import numpy as np
import scipy.sparse
from reac_diff_solver.conjgrad import conjugate_gradients


class Solver:
    def __init__(self, xBounds, yBounds, gridSize, initial_condition_function):
        """
        Initialises a Solver with the given bounds, gridsize and initial conditions. By default it has no reaction
        terms. These can be set using set_reactionFunction
        :param xBounds: x-range of the problem
        :type xBounds: list of two floats
        :param yBounds: y-range of the problem
        :type yBounds: list of two floats
        :param gridSize: number of gridpoints to use (in both the x- & y-direction)
        :type gridSize: int
        :param initial_condition_function: calculates the values of u and v at t=0 at each gridpoint
        :type initial_condition_function: function that takes two 2d numpy arrays and returns a list of two 2d numpy arrays
        """

        self.initial_condition_function = initial_condition_function
        self.set_grid(xBounds,yBounds,gridSize)                                         # create grid
        self.initialConditions_u = self.initial_condition_function(self.X, self.Y)[0]   # calculate initial conditions at grid points
        self.initialConditions_v = self.initial_condition_function(self.X, self.Y)[1]

        self.reactionFunction = lambda u, v, parameters: [0,0]                          # default reaction function returns 0 everywhere
        self.timeStepLength = 0.1                                                       # default time step

    def set_grid(self, xBounds, yBounds, gridSize):
        """
        Sets the grid parameters.
        :param xBounds: x-range of the problem
        :type xBounds: list of two floats
        :param yBounds: y-range of the problem
        :type yBounds: list of two floats
        :param gridSize: number of gridpoints to use (in both the x- & y-direction)
        :type gridSize: int
        :return:
        """
        self.xBounds = xBounds
        self.yBounds = yBounds
        self.gridSize = gridSize
        self.xStepLength = (xBounds[1] - xBounds[0])/(gridSize)                 # calculate x step length
        self.yStepLength = (yBounds[1] - yBounds[0])/(gridSize)                 # calculate y step length
        self.x = np.linspace(self.xBounds[0], self.xBounds[1], gridSize+1)[:-1] # generate coordinates of grid points & get rid of one of the boundary points
        self.y = np.linspace(self.yBounds[0], self.yBounds[1], gridSize+1)[:-1] # generate coordinates of grid points & get rid of one of the boundary points
        self.X, self.Y = np.meshgrid(self.x, self.y)                            # create mesh
        self.initialConditions = self.initial_condition_function(self.X, self.Y)# calculate initial conditions on new mesh

    def set_reactionFunction(self, function):
        """
        Set the reaction term of the equation (defaults to zero).
        :param function: calculates the value of the reaction terms at the given u and v
        :type function: function that takes two numpy arrays (containing values of u and v) and a list of parameters and returns a list of two numpy arrays
        :return:
        """
        self.reactionFunction = function

    def set_initialConditions(self, initial_condition_function):
        """
        Set the initial conditions used by the solver.
        :param initial_condition_function: calculates the values of u and v at t=0 at each gridpoint
        :type initial_condition_function: function that takes two 2d numpy arrays and returns a list of two 2d numpy arrays
        :return:
        """
        self.initial_condition_function = initial_condition_function
        self.initialConditions_u = self.initial_condition_function(self.X, self.Y)[0]   # calculate value of the initial conditions for u at the grid points
        self.initialConditions_v = self.initial_condition_function(self.X, self.Y)[1]   # calculate value of the initial conditions for v at the grid points

    def set_timeStepLength(self, length):
        """
        Set the size of the timestep used by the solver.
        :param length: timestep
        :type length: float
        :return:
        """
        self.timeStepLength = length

    def _create_arrays(self, times):
        """
        Create arrays to store solution at given times
        :param times: times at which to find the solution
        :type times: list of floats
        """

        self.uSolution = np.zeros((len(times), self.gridSize, self.gridSize))
        self.vSolution = np.zeros((len(times), self.gridSize, self.gridSize))


    def _create_fdmatrix(self):
        """
        Creates the finite difference matrix for periodic boundary conditions.
        """
        n = self.gridSize
        e = np.ones(n)
        diagonals = [e, -2 * e, e]
        offsets = [-1, 0, 1]
        L = scipy.sparse.spdiags(diagonals, offsets, n, n, format='csr').tolil()    # 1d laplacian with no boundary conditions

        # periodic boundary conditions
        L[0, -1] = 1
        L[-1, 0] = 1

        I = scipy.sparse.eye(n)

        # calculate 2d laplacian using kronecker product
        self.laplacian = scipy.sparse.kron(I, self.xStepLength**-2 * L) + scipy.sparse.kron(self.yStepLength**-2 * L, I)


    def solve(self, times, parameters, printProgress = False):
        """
        Solves the equation at the given times.
        :param times: times at which the solution is desired.
        :type times: list of floats
        :param parameters: parameters to give to the reaction function
        :type parameters: list
        :param printProgress: whether or not to print progress updates
        :return: the solution of u and v at the given times
        :rtype: list of two 2d numpy arrays (self.gridSize by self.gridSize)
        """
        # first two parameters are diffusion coefficients
        self.diff_u = parameters[0]
        self.diff_v = parameters[1]

        self._create_arrays(times)
        self._create_fdmatrix()

        self.uSolution[0] = self.initialConditions_u    # set first entry in solution array to initial condition
        self.vSolution[0] = self.initialConditions_v

        # reshape the initial conditions to create the initial vector to start the integration
        uvec = self.uSolution[0].reshape(-1)
        vvec = self.vSolution[0].reshape(-1)

        I = scipy.sparse.eye(self.gridSize**2)
        matrix_u = (I - self.timeStepLength * self.diff_u * self.laplacian) # iteration matrix for u
        matrix_v = (I - self.timeStepLength * self.diff_v * self.laplacian) # iteration matrix for v

        t = times[0] # starting time
        for i in range(1,len(times)):
            j = 1   # tracks number of steps in the iteration

            if printProgress:
                print("t=" + str(times[i]))

            # integrate until the next time in the times array
            while t < times[i]:
                # solve the IMEX equations for u and v, use the values at the previous iteration as the starting point
                # for the conjugate gradients function, as they should be close to the new value
                uvec = conjugate_gradients(matrix_u, uvec + self.timeStepLength * self.reactionFunction(uvec, vvec, parameters[2:])[0], uvec)[0]
                vvec = conjugate_gradients(matrix_v, vvec + self.timeStepLength * self.reactionFunction(uvec, vvec, parameters[2:])[1], vvec)[0]

                # update t & j
                t = times[i-1] + j*self.timeStepLength
                j += 1

            # reshape solutions back into 2d arrays of the values at each gridpoint
            self.uSolution[i] = uvec.reshape(self.gridSize, self.gridSize)
            self.vSolution[i] = vvec.reshape(self.gridSize, self.gridSize)

        return [self.uSolution, self.vSolution]
