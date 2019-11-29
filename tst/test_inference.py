import unittest
import numpy as np
from scipy.sparse import lil_matrix, csc_matrix
from scipy.optimize import rosen
from reac_diff_solver.solver import Solver
from reac_diff_solver.inference import Inference

class InferenceTest(unittest.TestCase):
    def test_set_model(self):
        testudata = np.zeros((2, 50, 50)) # doesn't matter
        testvdata = [0] # doesn't matter
        testtimes = [0] # doesn't matter
        infr = Inference(testudata, testvdata, testtimes)
        s = Solver(np.random.random(2), np.random.random(2), 50, lambda X, Y: [X, Y])
        infr.set_model(s.xBounds, s.yBounds, s.initial_condition_function)
        #self.assertTrue(np.array_equal(s.xBounds, infr.solver.xBounds))
        #self.assertTrue(np.array_equal(s.yBounds, infr.solver.yBounds))
        #self.assertEqual(s.gridSize, infr.solver.gridSize)
        s_vars = vars(s)
        infr_vars = vars(infr.solver)
        self.assertTrue(s_vars.keys() == infr_vars.keys())
        # Compare the two objects attribute by attribute
        for key in ["initial_condition_function", "initialConditions_u", "initialConditions_v", "xBounds", "yBounds", "gridSize", "xStepLength", "yStepLength", "x", "y", "X", "Y", "initialConditions"]:
            if isinstance(s_vars[key], (list, tuple, np.ndarray)):
                self.assertTrue(np.array_equal(s_vars[key], infr_vars[key]))
            else:
                self.assertTrue(s_vars[key] == infr_vars[key])

    def test_set_reaction_function(self):
        testudata = np.zeros((2, 50, 50)) # doesn't matter
        testvdata = [0] # doesn't matter
        testtimes = [0] # doesn't matter
        testfun = lambda X: X*X
        infr = Inference(testudata, testvdata, testtimes)
        s = Solver(np.random.random(2), np.random.random(2), 50, lambda X, Y: [X, Y])
        infr.solver = s
        infr.set_reaction_function(testfun)
        self.assertEqual(infr.solver.reactionFunction, testfun)
    
    def test_cost_func(self):
        N = 10
        testudata = np.array([np.tile(np.arange(1,N+1,1),[N, 1]), np.tile(np.arange(1,N+1,1),[N, 1])])
        testvdata = np.array([np.tile(np.arange(1,N+1,1),[N, 1]), np.tile(np.arange(1,N+1,1),[N, 1])])
        testtimes = np.array([0.0, 1.0])
        testparams = np.array([0.01, 0.03])
        infr = Inference(testudata, testvdata, testtimes)
        infr.set_model([0.0, 1.0], [0.0, 1.0], lambda X, Y: [X, Y])
        result = infr.cost_func(np.array([0.01, 0.03]))
        self.assertAlmostEqual(result, 13273.60382091)
    
    def test_fit_params(self):
        testudata = np.zeros((2, 50, 50)) # doesn't matter
        testvdata = [0] # doesn't matter
        testtimes = [0] # doesn't matter
        testx0 = [1.1, 1.2, 1.3, 1.4, 1.5]
        infr = Inference(testudata, testvdata, testtimes)
        infr.cost_func = rosen
        result = infr.fit_params(testx0)
        self.assertTrue(np.allclose(result, np.array([1.0, 1.0, 1.0, 1.0, 1.0]), rtol=1e-03))

if __name__ == "__main__":
    unittest.main()