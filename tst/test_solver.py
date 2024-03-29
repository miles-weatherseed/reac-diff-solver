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


    def test_solve_exponential(self):
        t = np.linspace(0, 5.0, 6)
        testparams = np.array([1.0, 2.0])
        k1, k2 = testparams

        def reaction_function(u, v, K):
            return [-K[0] * u, - K[1] * v]

        def exact(x, y, t):
            return [1.0 * np.exp(-k1 * t), 2.0 * np.exp(-k2 * t)]

        def initial_conditions(x, y):
            return exact(x, y, 0)

        solver = Solver([0.0, 1.0], [0.0, 1.0], 30, initial_conditions)
        solver.set_timeStepLength(0.0001)
        solver.set_reactionFunction(reaction_function)
        u, v = solver.solve(t, [5.0, 1.0, k1, k2])
        T, Y, X = np.meshgrid(t, solver.y, solver.x, indexing='ij')
        solution = exact(X, Y, T)

        error_u = np.max(np.abs((u - solution[0])))
        error_v = np.max(np.abs((v - solution[1])))
        print(error_u)
        print(error_v)

        # Check Uniform
        for i in range(len(t)):
            self.assertAlmostEqual(np.std(u[i]), 0)
            self.assertAlmostEqual(np.std(v[i]), 0)



        self.assertLess(error_u, 1e-4)
        self.assertLess(error_v, 1e-4)

    def test_solve_1ddiffusion(self):
        mode = 1.0
        D1 = 0.05
        D2 = 0.0
        t = np.linspace(0, 1.0, 11)

        def exact(x, y, t):
            return [np.sin(4 * np.pi * mode * x) * np.exp(-D1 * t * (4 * np.pi * mode) ** 2), 1]

        def initial_conditions(x, y):
            return exact(x, y, 0)

        solver = Solver([0.0, 1.0], [0.0, 1.0], 200, initial_conditions)
        solver.set_timeStepLength(0.0001)

        u, v = solver.solve(t, [D1, D2, 0])

        T, Y, X = np.meshgrid(t, solver.y, solver.x, indexing='ij')
        solution = exact(X, Y, T)


        error_u = np.max(np.abs((u - solution[0])))
        error_v = np.max(np.abs((v - solution[1])))

        # Check Uniform
        for i in range(len(t)):
            for j in range(len(solver.x)):
                self.assertAlmostEqual(np.std(u[i,:,j]), 0)

        self.assertAlmostEqual(np.std(v), 0)
        self.assertAlmostEqual(v[1,1,1], 1)

        self.assertLess(error_u, 3*1e-4)
        self.assertLess(error_v, 1e-4)


if __name__ == "__main__":
    unittest.main()