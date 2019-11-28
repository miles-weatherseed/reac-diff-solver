import numpy as np
import numpy.linalg as la

def conjugate_gradient(A, b, x0, tol = 0.0001):
    """
    Uses the conjugate gradient algorithm to solve the system Ax = b.
    :param A: matrix
    :param b: right-hand side
    :param x0: starting point
    :param tol: tolerance, iteration terminates when the norm of the residual has been reduced by the tolerance
    :return: solution of the equation, list of the norm of the residual after each iteration
    """

    x = np.copy(x0)     # solution
    r = b - A@x         # residual
    res = [la.norm(r)]  # stores the norm of the residuals at each step
    p = r               # initial search direction

    while res[-1] > tol * res[0]:
        alpha = np.dot(r,r)/np.dot(p,A@p)   # step length
        x += alpha * p                      # update x
        r_new = r - alpha * A@p             # new residual, store in new variable, because old r is needed to calculate beta
        beta = np.dot(r_new,r_new)/np.dot(r,r)
        r = r_new                           # update residual
        res.append(la.norm(r))
        p = r + beta * p                    # new search direction

    return x, res
