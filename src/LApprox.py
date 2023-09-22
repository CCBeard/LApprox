import numpy as np
import pandas as pd
import radvel
from scipy import spatial



def NDeriv_2(func, x0, dim1, dim2, **kwargs):
    '''
    Numerically calculate the second derivative of a function
    func: python function
    x0: point at which to calculate the derivative, N dimensional
    dim1: which dimension to take the first partial derivative of
    dim2: which dimension to take the second partial derivative of
    optional_args: list of optional arguments for the function
    eps: epsilon, the width of the numerical derivative estimate. Make sure it's smaller than all the parameters
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
    '''
    Function to calculate the Hessian Matrix of an RV likelihood
    '''
    Hessian = np.zeros([len(vals), len(vals)])
    for i in range(len(vals)):
        for j in range(len(vals)):
            Hessian[i][j] += (NDeriv_2(func, vals, i, j, **kwargs))
    print(Hessian)

    return np.array(Hessian)



def Laplace_Approximation(func, x0, **kwargs):
    '''
    Calculate the Laplace Approximation:

        Z = Integral(Exp(f(x))dx)

    Can be estimate as approximately:

        [(2pi)^2/(|det(H(x0)))]^(1/2) * exp(f(x0))

    Arguments:

        func: the function, f(x), in the exponent of the term we wish to estimate

        x0: the local maximum around which to compute the approxmation

    Returns:

        A: exp(f(x0))

        B: The term added by Laplace's approximation involving the Hessian Matrix


    '''

    logA = func(x0, **kwargs)


    H = Calculate_Hessian(func, x0, **kwargs)
    logB = np.log(((np.pi*2)**2/(np.abs(np.linalg.det(H))))**(1/2))
    print('Determinant: {}'.format(np.linalg.det(H)))



    return logA, logB
