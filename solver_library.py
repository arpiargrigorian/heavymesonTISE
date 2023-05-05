'''
Library of Solving Methods:

Using the Numerov Method, these solvers find the forwards/backwards solutions for a variety of potentials V(r) for the TISE: (d^2/dr^2 - V(r))u(r) = Eu(r). 
The solution u(r) is actually the reduced solution from the radial wavefunction R(r). u(r)/r = R(r). 

This library needs numerov_library.py in the same directory to work. 

'''
import numpy as np #Established Modules
from scipy.special import airy

from numerov_library import * #My Module

def centrf(rho, l):
    '''
    Calculates the value of the centrifugal potential l(l + 1)/rho^2
    
    Inputs:
      rho is a number
      l is a number; it should be a positive integer to give a physically relevant answer
      
     Returns:
      The value of the centrifugal potential
    
    '''
    return l*(l + 1)/rho**2

def Coulomb(rho, coulomb_constant):
    '''
    Calculates the coulomb potential coulomb_constant/rho
    
    Inputs:
      rho is a number
      coulomb constant is a positive number
     Returns:
      The value of the coulomb potential
    '''
    return coulomb_constant/rho

def numerov_backward_centrf_lin_solver_v2(el, coulomb_constant, rho, zeta, guess):
 
    '''
    Finds the backwards moving solution to the equation:
      (d^2/drho^2 - l(l + 1)/rho^2 - rho + coulomb_constant/rho + \zeta)u(r) = 0
    
    Inputs:
      el: a positive integer. It is the angular momentum quantum number.
      coulomb_constant: A positive number. It is the strength of the frozen Coulomb potential
      rho: a numpy array of ordered numbers. It should correspond to the region where the numerov method will not explode. 
      zeta: a number that corresponds to the guess for the energy eigen value. 
      guess: The last three values of the solution, given as an ordered numpy array or list. This is the initial condition for the problem. 
    
    Returns:
      h: an array containting the spacing between h[i] = \Delta rho[i] = rho[i + 1] - rho[i]
      backwards_sol: an array of numbers of length len(rho) corresponding to the solution to the equation. 
      
    ''' 
    
    ## Initialize Arrays
    rhonum = len(rho)
    backward_sol = np.zeros(rhonum)
    
    ## Assign the first two initial conditions to the array
    backward_sol[-1] = guess[-1]
    backward_sol[-2] = guess[-2]
    backward_sol[-3] = guess[-3]
    
    h = rho[1:] - rho[:-1]
    
    ## Solve:
    
    for i in range(3, rhonum):
        
        F = np.zeros(3)
        Y = np.zeros(2)
        
        F[0] = -1*(rho[-1*i - 1] - zeta) - centrf(rho[-1*i - 1], el) + Coulomb(rho[-1*i - 1], coulomb_constant)
        F[1] = -1*(rho[-1*i] - zeta)- centrf(rho[-1*i], el) + Coulomb(rho[-1*i], coulomb_constant)
        F[2] = -1*(rho[-1*i + 1] - zeta)- centrf(rho[-1*i + 1], el) + Coulomb(rho[-1*i + 1], coulomb_constant)
        
        Y[0] = backward_sol[-1*i]
        Y[1] = backward_sol[-1*i + 1]
        
        backward_sol[-1*i - 1] = numerov_backward_step(F, Y, h[-1*i])
        
    return h, backward_sol


def backward_centrf_linear_matrix(el, coulomb_constant, rho, zeta_array):
    '''
    Creates a matrix where each row is the backwards moving solution for the TISE with the Cornell Potential and zeta in zeta_array
  
    Inputs:
      el: a positive integer. It is the angular momentum quantum number.
      coulomb_constant: A positive number. It is the strength of the frozen Coulomb potential
      rho: a numpy array of ordered numbers. It should correspond to the region where the numerov method will not explode. 
      zeta_array: an ordered array of numbers that corresponds to a range of zeta values. 
    
    Returns:
      backwards_solutions: a Matrix where each row i corresponds to the backwards moving solution over rho for zeta value zeta_array[i] 
      
    '''
    backward_solutions = np.zeros((len(zeta_array), len(rho)))
    for i, zeta in enumerate(zeta_array):
        Ai, A, B, C = airy(rho - zeta)
        hB, backward_solutions[i] = numerov_backward_centrf_lin_solver_v2(el, coulomb_constant, rho, zeta, [Ai[-3], Ai[-2], Ai[-1]])
    return backward_solutions
  
def numerov_forward_centrf_lin_solver_v2(el, coulomb_constant, rho, zeta, guess):
 
    '''
    Finds the forwards moving solution to the equation:
      (d^2/drho^2 - l(l + 1)/rho^2 - rho + coulomb_constant/rho + \zeta)u(r) = 0
    
    Inputs:
      el: a positive integer. It is the angular momentum quantum number.
      coulomb_constant: A positive number. It is the strength of the frozen Coulomb potential.
      rho: a numpy array of ordered numbers. It should correspond to the region where the numerov method will not explode. 
      zeta: a number that corresponds to the guess for the energy eigen value. 
      guess: The first three values of the solution, given as an ordered numpy array or list. This is the initial condition for the problem. 
    
    Returns:
      h: an array containting the spacing between h[i] = \Delta rho[i] = rho[i + 1] - rho[i]
      forward_sol: an array of numbers of length len(rho) corresponding to the forward (right) moving solution to the equation. 
    ''' 
    
    ## Initialize Arrays
    rhonum = len(rho)
    forward_sol = np.zeros(rhonum)
    
    ## Assign the first two initial conditions to the array
    forward_sol[0] = guess[0]
    forward_sol[1] = guess[1]
    forward_sol[2] = guess[2]
    
    h = rho[1:] - rho[:-1]
    
    ## Solve:
    
    for i in range(2, rhonum - 1):
        
        F = np.zeros(3)
        Y = np.zeros(2)
        
        F[0] = -1*(rho[i - 1] - zeta) - centrf(rho[i - 1], el) + Coulomb(rho[i - 1], coulomb_constant)
        F[1] = -1*(rho[i] - zeta) - centrf(rho[i], el) + Coulomb(rho[i], coulomb_constant)
        F[2] = -1*(rho[i + 1] - zeta) - centrf(rho[i + 1], el) + Coulomb(rho[i + 1], coulomb_constant)
        
        Y[0] = forward_sol[i - 1]
        Y[1] = forward_sol[i]

        
        forward_sol[i + 1] = numerov_forward_step(F, Y, h[i])
        
    return h, forward_sol

def forward_centrf_linear_matrix(el, coulomb_constant, rho, zeta_array, n = 1):
    '''
    Creates a matrix where each row is the forwards moving solution to the TISE with a Cornell potential for zeta in zeta_array
  
    Inputs:
      el: a positive integer. It is the angular momentum quantum number.
      coulomb_constant: A positive number. It is the strength of the frozen Coulomb potential
      rho: a numpy array of ordered numbers. It should correspond to the region where the numerov method will not explode. 
      zeta_array: an ordered array of numbers that corresponds to a range of zeta values. 
      n: the principal quantum number. Should be a positive integer >= 1. 
    
    Returns:
      forward_solutions: a Matrix where each row i corresponds to the forwards moving solution over rho for zeta value zeta_array[i] 
    ''' 
    forward_solutions = np.zeros((len(zeta_array), len(rho)))
    for i, zeta in enumerate(zeta_array):
        sign = (-1)**(n + 1)
        u_0 = [sign*rho[0]**(el + 1), sign*rho[1]**(el + 1), sign*rho[2]**(el + 1)]
        hF, forward_solutions[i] = numerov_forward_centrf_lin_solver_v2(el, coulomb_constant, rho, zeta, u_0)
    return forward_solutions
