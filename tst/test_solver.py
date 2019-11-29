import unittest
import numpy as np
from reac_diff_solver.solver import Solver


class SolverTest(unittest.TestCase):
    def test_set_grid(self):
        testxbounds = np.random.random(2)
        testybounds = np.random.random(2)
        testgridsize = 50
        testfun = lambda X,Y: [np.sin(4*np.pi*X) * np.sin(2*np.pi*Y), np.sin(4*np.pi*X) * np.sin(2*np.pi*Y)]
        # Expected values
        expected_x = np.linspace(testxbounds[0], testxbounds[1], testgridsize + 1)[:-1]
        expected_y = np.linspace(testybounds[0], testybounds[1], testgridsize + 1)[:-1]
        expected_X, expected_Y = np.meshgrid(expected_x, expected_y)
        # Test
        s = Solver(np.random.random(2), np.random.random(2), np.random.randint(1,100), testfun)
        s.set_grid(testxbounds, testybounds, testgridsize)
        # Checks
        self.assertTrue(np.array_equal(s.xBounds, testxbounds))
        self.assertTrue(np.array_equal(s.yBounds, testybounds))
        self.assertEqual(s.gridSize, testgridsize)
        self.assertEqual(s.xStepLength, (testxbounds[1] - testxbounds[0])/(testgridsize))
        self.assertEqual(s.yStepLength, (testybounds[1] - testybounds[0])/(testgridsize))
        self.assertTrue(np.array_equal(s.x, expected_x))
        self.assertTrue(np.array_equal(s.y, expected_y))
        self.assertTrue(np.array_equal(s.X, expected_X))
        self.assertTrue(np.array_equal(s.Y, expected_Y))
        self.assertTrue(np.array_equal(s.initialConditions, testfun(expected_X, expected_Y)))

    def test_set_reactionFunction(self):
        testfun = lambda X: 2*X
        s = Solver(np.random.random(2), np.random.random(2), np.random.randint(1,100), lambda X, Y: [X, Y])
        s.set_reactionFunction(testfun)
        self.assertEqual(s.reactionFunction, testfun)

    def test_set_initialConditions(self):
        testfun = lambda X,Y: [np.sin(X) * np.sin(Y), np.cos(X) * np.cos(Y)]
        s = Solver(np.random.random(2), np.random.random(2), np.random.randint(1,100), lambda X, Y: [X, Y])
        s.set_initialConditions(testfun)
        self.assertTrue(np.array_equal(s.initialConditions_u, testfun(s.X, s.Y)[0]))
        self.assertTrue(np.array_equal(s.initialConditions_v, testfun(s.X, s.Y)[1]))
        self.assertEqual(s.initial_condition_function, testfun)

    def test_set_timeStepLength(self):
        testlen = np.random.random(1)
        s = Solver(np.random.random(2), np.random.random(2), np.random.randint(1,100), lambda X, Y: [X, Y])
        s.set_timeStepLength(testlen)
        self.assertAlmostEqual(s.timeStepLength, testlen)

    def test_solve(self):
        testtimes = np.array([0.0, 1.0])
        testparams = np.array([0.01, 0.03])
        s = Solver([0.0, 1.0], [0.0, 1.0], 100, lambda X, Y: [X, Y])
        u, v = s.solve(testtimes, testparams)
        self.assertTupleEqual(u.shape, (2, 100, 100))
        self.assertTupleEqual(v.shape, (2, 100, 100))

    def test_solve_oscillation(self):

        def uniform_initial_conditions(X,Y):
            """
            Creates uniform initial conditions with u = 1, v = 0
            """
            return [np.ones_like(X),np.zeros_like(X)]

        def oscillation_reactionfunction(u,v, parameters):
            """
            Reaction function with a sinusoidal exact solution.
            """
            w = parameters[0] # angular speed
            return [w*v, -w*u]

        testtimes = np.linspace(0,2*np.pi, 20)
        for w in [1.0, 2.0, 3.0]: # angular speed
            testparams = np.array([0.01, 0.01, w])
            s = Solver([0.0, 1.0], [0.0, 1.0], 20, uniform_initial_conditions)
            s.set_reactionFunction(oscillation_reactionfunction)
            s.set_timeStepLength(0.0001)
            u, v = s.solve(testtimes, testparams)

            for i in range(len(testtimes)):
                # solution remains uniform
                self.assertAlmostEqual(np.std(u[i]), 0)
                self.assertAlmostEqual(np.std(v[i]), 0)

                # solution is sine & cosine with given angular speed
                self.assertAlmostEqual(u[i,10,10], np.cos(w * testtimes[i]), places = 3)
                self.assertAlmostEqual(v[i,10,10], -np.sin(w * testtimes[i]), places = 3)

    def test_solve_2dheatequation(self):
        def sinusoidal_initial_conditions(X,Y):
            """
            Creates initial conditions with a sinusoidal wave in the x and y directions.
            """
            return [np.sin(4*np.pi*X) * np.sin(2*np.pi*Y), np.sin(4*np.pi*X) * np.sin(2*np.pi*Y)]

        testtimes = np.linspace(0.0, 1.0, 20)
        testparams = [0.05, 0.01] # diffusion coefficients
        s = Solver([0.0, 1.0],[0.0, 1.0], 50, sinusoidal_initial_conditions)
        s.set_timeStepLength(0.0001)
        u, v = s.solve(testtimes, testparams)

        for i in range(len(testtimes)):
            for j in range(50):
                for k in range(50):
                    self.assertAlmostEqual(u[i,j,k], np.sin(4*np.pi*s.X[j,k]) * np.sin(2*np.pi*s.Y[j,k]) * np.exp(-testparams[0]*20*np.pi**2 * testtimes[i]), places = 2)
                    self.assertAlmostEqual(v[i,j,k], np.sin(4*np.pi*s.X[j,k]) * np.sin(2*np.pi*s.Y[j,k]) * np.exp(-testparams[1]*20*np.pi**2 * testtimes[i]), places = 2)

if __name__ == "__main__":
    unittest.main()