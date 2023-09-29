import numpy as np
import pandas as pd
import radvel
from scipy import spatial



def NDeriv_2(func, x0, dim1, dim2, **kwargs):
    '''Numerically calculate the second partial derivative of a function

    Args:
        func (function): python function to calculate the derivitive of with dimension N
        x0 (array): point at which to calculate the derivative, N dimensional
        dim1 (int): which dimension to take the first partial derivative of
        dim2 (int): which dimension to take the second partial derivative of
        kwargs: keyword arguments

    Return:
        float: numerical second partial derivative with respect to dimensions one and two
    '''

    try:
        eps = kwargs['eps']
    except KeyError:
        eps = 1e-5


    x_upup = x0.copy()
    x_updown = x0.copy()
    x_downup = x0.copy()
    x_downdown = x0.copy()
    x_upsame = x0.copy()
    x_downsame = x0.copy()

    #scale the parameter adjustment by dimension, but this fails if it is exactly zero
    dim1adj = eps/2
    dim2adj = eps/2


    #For mixed second derivitives
    x_upup[dim1] += dim1adj
    x_upup[dim2] += dim2adj

    x_updown[dim1] += dim1adj
    x_updown[dim2] -= dim2adj

    x_downup[dim1] -= dim1adj
    x_downup[dim2] += dim2adj

    x_downdown[dim1] -= dim1adj
    x_downdown[dim2] -= dim2adj

    func_upup = func(x_upup, **kwargs)
    func_updown = func(x_updown, **kwargs)

    func_downup = func(x_downup, **kwargs)
    func_downdown = func(x_downdown, **kwargs)

    #For un-mixed second derivitives
    x_upsame[dim1] += dim1adj  #probably unecessary
    x_downsame[dim1] -= dim1adj  #probably unecessary

    func_upsame = func(x_upsame, **kwargs) #probably unecessary
    func_downsame = func(x_downsame, **kwargs) #probably unecessary
    func_x0 = func(x0, **kwargs) #probably unecessary

    #if dim1 == dim2:
        #derivitive is not mixed, both are in same dimension
        #out = (2*func_x0 - func_upsame - func_downsame)/(dim1adj)**2
    #else:
    out = (func_upup - func_updown - func_downup + func_downdown)/(4*dim1adj*dim2adj)

    return out

def Calculate_Hessian(func, vals, **kwargs):
    '''Function to calculate the Hessian Matrix of a generic function

    Args:
        func (function): Python function to estimate the Hessian matrix of, N dimensional
        vals (array): numpy array of length N. The point at which to estimate the matrix of second derivatives
        kwargs: keyword arguments
    Returns:
        array: numpy array, Hessian of the function func. Shape NxN.
    '''
    Hessian = np.zeros([len(vals), len(vals)])
    for i in range(len(vals)):
        for j in range(len(vals)):
            Hessian[i][j] += (NDeriv_2(func, vals, i, j, **kwargs))

    return np.array(Hessian)



def Laplace_Approximation(func, x0, **kwargs):
    r'''Calculate the Laplace Approximation of an integral of specific form.

    A challenging integral, when possible to write in terms of an exponent:
        .. math::
            Z = \int(exp(f(x))dx)

    Can be estimated as approximately:

        .. math::
            [\frac{(2\pi)^{2}}{det|H(x_{0})|}]^{\frac{1}{2}} * exp(f(x_{0}))

    where H is the functions Hessian matrix, and x0 is a region of high probability.

    Arguments:
        func (function): python function, f(x), in the exponent of the term we wish to estimate. Note that this is NOT
            the total function that we are trying to integrate, but f(x), in the exponent.

        x0 (array): numpy array of N values, where the function is N-dimensional. This is the
            local maximum around which to compute the approxmation. When done correctly, one should optimize a function
            first before finding the Laplace Approximation.

    Returns:
        (tuple): logA (float) and logB (float). logA = f(x0), where exp(f(x0)) is A. We return the logarithm because, in practice, many functions we wish
        to calculate the Laplace Approximation for return very small or very large values that can overflow a computer's
        floating point precision. The term added by Laplace's approximation involving the Hessian Matrix, [(2\pi)^2/(det|H(x_{0})|)]^1/2
        where we take the logarithm for consistency with A.
    '''

    logA = func(x0, **kwargs)


    H = Calculate_Hessian(func, x0, **kwargs)
    logB = np.log(((np.pi*2)**2/(np.abs(np.linalg.det(H))))**(1/2))
    print('Determinant: {}'.format(np.linalg.det(H)))



    return logA, logB
