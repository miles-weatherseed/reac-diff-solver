import unittest
import numpy as np
from solver import Solver

class TestModel(unittest.TestCase):
    xbounds = [0.0, 2.0]
    ybounds = [0.0, 4.0]
    gridsize = 4

    def initial_function(self, X, Y):
        return np.sin(4 * np.pi * X) * np.sin(2 * np.pi * Y), np.sin(4 * np.pi * X) * np.sin(2 * np.pi * Y) #* np.exp(-40 * np.pi ** 2 * T)

    def exact_solution(self, X, Y, T):
        return np.sin(4 * np.pi * X) * np.sin(2 * np.pi * Y) * np.exp(-40 * np.pi ** 2 * T),

    def u_reaction(self, u, v):
        return 0

    def v_reaction(self, u, v):
        pass

    def setUp(self) -> None:
        self.sol = Solver(self.xbounds, self.ybounds, self.gridsize, self.initial_function)

    def test_init(self):
        self.assertEqual(self.sol.xBounds, self.xbounds)
        self.assertEqual(self.sol.yBounds, self.ybounds)
        self.assertEqual(self.sol.gridSize, self.gridsize)
        self.assertEqual(self.sol.xStepLength, 0.5)
        self.assertEqual(self.sol.yStepLength, 1.0)
        self.assertEqual(self.sol.x[1]- self.sol.x[0], self.sol.xStepLength)
        self.assertIsNone(np.testing.assert_array_equal(self.sol.x, np.arange(self.xbounds[0], self.xbounds[1], 0.5)))
        self.assertIsNone(np.testing.assert_array_equal(self.sol.y, np.arange(self.ybounds[0], self.ybounds[1], 1.0)))
        print(self.sol.x)

    def test_matrix(self):
        self.sol.solve([0.1, 0.2, 0.4], [.5, 1.0])
        print(self.sol.laplacian)

    def test_set_bounds(self):
        xbounds = [10, 20]
        self.sol.set_xBounds(xbounds)
        self.assertEqual(self.sol.xBounds, xbounds)

        ybounds = [10, 20]
        self.sol.set_yBounds(ybounds)
        self.assertEqual(self.sol.yBounds, ybounds)



