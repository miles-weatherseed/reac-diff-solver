import unittest
import numpy as np
from scipy.sparse import lil_matrix, csc_matrix
from reac_diff_solver.conjgrad import conjugate_gradients

class ConjugateGradientsTest(unittest.TestCase):
    def test_non_square(self):
        A = np.arange(0,6,1).reshape(3,2)
        b = np.zeros(3)  # Don't care about validity
        x0 = np.zeros(3) # Don't care about validity
        with self.assertRaises(TypeError) as arctx:
            conjugate_gradients(A, b, x0)
        self.assertEqual(str(arctx.exception), "The matrix provided is not square.")

    def test_dim_mismatch_b(self):
        N = 3
        A = np.arange(0,N*N,1).reshape(N,N)
        b = np.zeros(N + 1)
        x0 = np.zeros(N) # Don't care about validity
        with self.assertRaises(TypeError) as arctx:
            conjugate_gradients(A, b, x0)
        self.assertEqual(str(arctx.exception), "The dimensions of the right hand side do not match the dimensions of the matrix provided.")

    def test_dim_mismatch_x0(self):
        N = 3
        A = np.arange(0,N*N,1).reshape(N,N)
        b = np.zeros(N)
        x0 = np.zeros(N + 1)
        with self.assertRaises(TypeError) as arctx:
            conjugate_gradients(A, b, x0)
        self.assertEqual(str(arctx.exception), "The dimensions of the starting point do not match the dimensions of the matrix provided.")
    
    def test_result_converging_dense(self):
        N = 3
        A = np.array([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])
        b = np.ones(N)
        x0 = np.zeros(N)
        results = conjugate_gradients(A, b, x0)
        self.assertTrue(np.allclose(results[0], np.array([1.5, 2.0, 1.5])))

    def test_result_converging_sparse(self):
        N = 1000
        A = lil_matrix((N,N))
        A.setdiag(np.arange(1,N+1,1))
        A = np.eye(N) + A*A.transpose()
        b = np.ones(N)
        x0 = np.zeros(N)
        results = conjugate_gradients(csc_matrix(A), b, x0, nmax = 2*N)
        expected_results = np.ones(N) / A.diagonal()
        self.assertTrue(np.allclose(results[0], expected_results))

    def test_result_nonconverging(self):
        N = 3
        A = np.array([[3, -1, 4], [5, 1, 8], [1, 2, 0]]) # not definite because one eigenvalue is negative
        b = np.ones(N)
        x0 = np.zeros(N)
        with self.assertRaises(Exception) as arctx:
            results = conjugate_gradients(A, b, x0)
        self.assertEqual(str(arctx.exception), "The iteration has failed to converge within nmax(=100) iterations.")
        # If solved the examples[0] would have been [1, 0, -0.5]

if __name__ == "__main__":
    unittest.main()