'''
The Numerov Library:
  Contains all functions needed to use the Numerov Method to Solve a 2nd Order ODE
'''

import numpy as np

## The Integration Function (equivalent to numpy.trapz(y, x=x))
def trap_v2(x, y):
    '''
    Integrate along the axis of x using the trapezoidal method. 
    
    Inputs:
      x is a 1D numpy array. x should be ordered.  
      y is a 1D numpy array representing the function y = f(x). The lengths of x and y must be identical. 
    
    Returns:
      total_area: the area under the curve from x[0] to x[-1]
    '''
    
    ### Exceptions
    if type(x) != np.ndarray:
      raise Exception("x must be a numpy array or list.")

    if type(y) != np.ndarray or type(y) != list:
      raise Exception("y must be a numpy array or list.")      
    
    if y.shape != x.shape:
      raise Exception("x and y do not have the same dimensions {} and {}".format(x.shape, y.shape))
    
    # Code
    h = x[1:] - x[:-1]
    area = h*(y[:-1] + y[1:])/2.
    total_area = np.sum(area)
    return total_area

### Numerov Method:
### Define all the terms in the full equation. 

def N_0(Fx, Yx, h):
    #The first term in the numerator for the numerov method.
    # Fx, Yx, and h are all numbers
    return 2*(1 - (5./12.)*h**2*Fx)*Yx

def N_D(Fx,h): 
    # This can be the denominator of the numerov method
    # Or the second term in the numerator
    # Fx, Yx, and h are all numbers
    return 1 + (1./12.)*h**2*Fx

def numerov_forward_step(F, Y, h):
    '''
    Take a forward (right moving) step to find the value of the solution Y(x + h) using the Numerov Method. 
    
    Inputs:
      F is a list or array that has exactly three elements (numbers) and corespond to the assigment F[0] = F(x - h), F[1] = F(x), F[2] = F(x + h)
      Y is a list or array with exactly two element (numbers), Y[0] =  Y(x - h), Y[1] = Y(x)
      h is a number that represents the step size
      
    Returns the value of Y(x + h) 
    '''
    return (N_0(F[1], Y[1], h) - N_D(F[0], h)*Y[0])/N_D(F[2], h)

def numerov_backward_step(F, Y, h):
    '''
    Take a backward (left moving) step to find the value of the solution Y(x - h) using the Numerov Method. 
    
    Inputs:
      F is a list or array that has exactly three elements (numbers) and corespond to the assigment F[0] = F(x - h), F[1] = F(x), F[2] = F(x + h)
      Y is a list or array with exactly two element (numbers), Y[0] =  Y(x), Y[1] = Y(x + h)
      h is a number that represents the step size
      
    Returns the value of Y(x - h) 
    '''
    return (N_0(F[1], Y[0], h) - N_D(F[2], h)*Y[1])/N_D(F[0], h)
