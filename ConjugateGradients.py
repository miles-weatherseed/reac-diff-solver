import numpy as np
import numpy.linalg as la

def conjugate_gradients(A, b, x0, tol = 0.0001, nmax = 100):
    """
    Uses the conjugate gradient algorithm to solve the system Ax = b.

    :param A: symmetric, positive definite matrix
    :type A: 2d numpy array or scipy sparse matrix
    :param b: right-hand side of the system
    :type b: 1d numpy array
    :param x0: starting point of the iteration
    :type x0: 1d numpy array
    :param tol: tolerance, iteration terminates when the norm of the residual has been reduced by the tolerance
    :type tol: float
    :param nmax: maximum number of iterations
    :type nmax: int
    :raises TypeError: the matrix A is not square
    :raises TypeError: the dimensions of A & b don't match
    :raises TypeError: the dimensions of A & x0 don't match
    :raises Exception: the iteration has failed to converge in nmax iterations
    :return: solution of the equation, list of the norm of the residual after each iteration
    :rtype: tuple of 1d numpy array and a list of floats
    """

    # check dimensions are compatible
    N = A.shape[0]
    if A.shape[1] != N:
        raise TypeError('The matrix provided is not square.')
    if b.shape[0] != N:
        raise TypeError('The dimensions of the right hand side do not match the dimensions of the matrix provided.')
    if x0.shape[0] != N:
        raise TypeError('The dimensions of the starting point do not match the dimensions of the matrix provided.')

    x = np.copy(x0)     # solution
    r = b - A@x         # residual
    res = [la.norm(r)]  # stores the norm of the residuals at each step
    p = r               # initial search direction

    n = 0               # number of iterations passed
    while res[-1] > tol * res[0]:
        alpha = np.dot(r,r)/np.dot(p,A@p)   # step length
        x += alpha * p                      # update x
        r_new = r - alpha * A@p             # new residual, store in new variable, because old r is needed to calculate beta
        beta = np.dot(r_new,r_new)/np.dot(r,r)
        r = r_new                           # update residual
        res.append(la.norm(r))
        p = r + beta * p                    # new search direction

        n+=1                                # update number of iterations
        if n > nmax:
            raise Exception('The iteration has failed to converge within nmax(=' + str(nmax) + ') iterations.')

    return x, res