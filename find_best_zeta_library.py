import timeit
from solver_library import *


def diff_ratio(rho, func):
    ## Finds func'/func
    func_prime = np.diff(func)/np.diff(rho)
    return func_prime/func[:-1]

def Deltasq_N(rho, forward_sol, backward_sol):
    '''
    Finds the Numerator of the \Delta^2 functional
    
    Input:
        rho: an ordered array of numbers.
        forward_sol: an array of numbers representing the forward moving solution u(rho)
        backward_sol: an array of numbers representing the backward moving solution u(rho)
        
    Returns the Numerator of \Delta^2
    
    '''
    return np.sum((diff_ratio(rho, forward_sol) - diff_ratio(rho, backward_sol))**2)

def Deltasq_D(rho, forward_sol, backward_sol):
    '''
    Finds the Denominator of the \Delta^2 functional
    
    Input:
        rho: an ordered array of numbers.
        forward_sol: an array of numbers representing the forward moving solution u(rho)
        backward_sol: an array of numbers representing the backward moving solution u(rho)
        
    Returns the Denomiator of \Delta^2    
    
    '''
    print(diff_ratio(rho, forward_sol))
    print(diff_ratio(rho, backward_sol))
    np.sum((diff_ratio(rho, forward_sol) + diff_ratio(rho, backward_sol))**2)

'''
Below, I will already have the rho_match cut to the right length and the matricies cut to the right length so it is
not something I need to test for. 

'''    

def sum_function_distance(forward_solution, backward_solution):
    '''
    Finds the squared difference between two 1D numpy arrays
    
    '''
    return np.sqrt(np.sum((forward_solution - backward_solution)**2))

def normalize_solution(rho, solution_matrix):
    
    '''
    Normalizes the solutions present in each row seperatley. 
    
    
    Inputs:
    
        rho: an ordered 1D numpy array. 
        
        solution_matrix: a 2D numpy array. Each row corresponds to the solution of the ODE with a different 
            eigenvalue \zeta. The columns represent the solution along the rho axis and the number of colummsn
            match the number of rho values. 
        
    
    Returns: 
    
        normalized_sol_matrix: The solution_matrix with each solution row properly normalized. 
    '''
    
    normalized_sol_matrix = np.zeros(solution_matrix.shape)
    
    for i in range(len(solution_matrix)):
        norm = np.sqrt(np.trapz(solution_matrix[i]**2, rho))
        normalized_sol_matrix[i] = solution_matrix[i]/norm
    
    return normalized_sol_matrix
        

def cut_backwards_sol_matrix(rho_index_max, sol_matrix):
    '''
    Slices the solution matrix to an upper index rho_index_max
    '''
    
    return sol_matrix[:, :rho_index_max]

def cut_forwards_sol_matrix(rho_index_min, sol_matrix):
    '''
    Slices the solution matrix to some lower index rho_index_min
    '''
    
    return sol_matrix[:, rho_index_min:]
    
def best_fit_zeta(rho_match, forward_solution_matrix_cut, backward_solution_matrix_cut):
    '''
    Finds which index corresponds to the solution with the best matching between forwards and backwards solutions. 
    
    Input:
        rho_match: an ordered numpy array. 
            a slice of the rho array that the forward and backward solutions have in common. 
        forward_solution_matrix_cut: A 2D numpy array. The row number corresponds to the number of zeta values 
            zeta_num, and the column number must be the same length as rho_match. This matrix corresponds to
            the forward moving solutions for different values of zeta.
        backward_solution_matrix: A 2D numpt matrix with the same shape as forwad_solution_matrix_cut. 
    Returns: 
        best_fit_index: the index corresponding to the row that minimized the Delta^2 function.
        best_fit_deltasq: The value of the minimum deltasq
        
    '''
    zeta_num = len(forward_solution_matrix_cut)
    delta_array = np.zeros(zeta_num)
    
    for i in range(zeta_num):
        
        forward_sol = forward_solution_matrix_cut[i]
        backward_sol = backward_solution_matrix_cut[i]
        
        delta_array[i] = sum_function_distance(forward_sol, backward_sol)

    best_fit_index = np.argmin(delta_array)
    best_fit_deltasq = np.min(delta_array)
    
    return best_fit_index, best_fit_deltasq

def zeta_interval_test(index, zeta_array):
    '''
    Tests whether the best fit index is on the edge of the zeta_array interval
    '''
    
    if index == 0:
        return False
    elif index == len(zeta_array) - 1: # On the edge
        return False
    else:
        return True

def accuracy(zeta_init_interval, zeta_num_array):
    zeta_min = zeta_init_interval[0]
    zeta_max = zeta_init_interval[1]
    for i in range(len(zeta_num_array)):
        zeta_space = np.linspace(zeta_min,zeta_max, num = zeta_num_array[i])
        h = np.diff(zeta_space)[0] # Just one of them
        print('The accuracy of iteration {} will be {} GeV'.format(i, h))
        zeta_min = zeta_space[0]
        zeta_max = zeta_space[2]
    return 0
    
def iter_find_zeta(el, quantum_n, coulomb_constant, rho, rho_cut_interval, zeta_init_interval, zeta_num_array, threshold = 10, SURPRESS_READOUT = False):
    '''
    Find the best fit zeta value in a range of zeta values for the 3D TISE Cornell Potential Problem
    
    Inputs:
        el: an integer >= 0. the total angular momentum quantum number.
        quantum_n: An interger >= 1. The principa Quantum Number
        coulomb_constant: The normalized strength of the Coulomb term. A positive real number. 
        rho: A 1D, ordered numpy array. the full array of rho values. 
        rho_cut_interval: a list or tuple of two integers that represent the minium value of the backwards solution
            and the maxiumum rho that the forwards moving solution reaches, respectively. 
        zeta_init_interval: A list or tuple of two numbers. The min/max guess for the zeta values. 
        zeta_num_array: A list or tuple of n positive integers. The length of the array is the number of times that the
            function will solve the TISE problem for an interval. The number is the number of zeta values per interval
            iteration. 
    Returns:
        zeta_value: The best value of zeta. 
        Deltasq: The value of delta sq that acted as the minimum. 
        best_backwards_sol: A 2XN matrix. best_backwards_sol[0] is the rho values for the best backwards solution.
            best_backwards_sol[1] is the solution u(rho).
        best_forwards_sol: A 2XN matrix. best_forwards_sol[0] is the rho values for the best forwards solution.
            best_forwards_sol[1] is the solution u(rho).
    '''
    
    ## Calculate the Accuracy of the zeta_num_array choice
    if SURPRESS_READOUT == False:
        throwaway =  accuracy(zeta_init_interval, zeta_num_array)
    
    rho_min_cut = rho_cut_interval[0]
    rho_max_cut = rho_cut_interval[1]
    
    rhoF = rho[:rho_max_cut] 
    rhoB = rho[rho_min_cut:]
    
    total_iter_num = len(zeta_num_array)
    
    # Set the first zeta_min and zeta_max
    zeta_min = zeta_init_interval[0]
    zeta_max = zeta_init_interval[1]
    
    start_time = timeit.default_timer()
    
    ## Initialize Some Arrays I will use
    best_fit_zeta_array = np.zeros(total_iter_num)
    best_fit_deltasq_array = np.zeros(total_iter_num)
    
    for i in range(total_iter_num):
        
        zeta_array = np.linspace(zeta_min, zeta_max, num = zeta_num_array[i])
        
        ## Solve the Equations:
        forward_solutions = forward_centrf_linear_matrix(el, coulomb_constant, rhoF, zeta_array, n = quantum_n)
        backward_solutions = backward_centrf_linear_matrix(el, coulomb_constant, rhoB, zeta_array)
        
        
        ## Normalize Each Solution:
        forward_normalized_solutions = normalize_solution(rhoF, forward_solutions)
        backward_normalized_solutions = normalize_solution(rhoB, backward_solutions)

        ## This did not work for some reason. 
        #for i in range(len(rhoF)):
        #    if any(sol > threshold for sol in forward_normalized_solutions[:, i]):
        #        raise Exception('The Forward Solution Exploded Past Threshold at rho = {}. Find other values of rho cut for the solution to converge'.format(rhoF[i]))
        #for i in range(len(rhoB)):
        #    if any(sol > threshold for sol in backward_normalized_solutions[:, i]):
        #        raise Exception('The Backward Solution Exploded Past Threshold at rho = {}. Find other values of rho cut for the solution to converge'.format(rhoB[i]))

        
        ## Cut The Matricies
        rho_match = rho[rho_min_cut:rho_max_cut]
        forward_match_matrix = cut_forwards_sol_matrix(rho_min_cut, forward_normalized_solutions )
        backward_match_matrix = cut_backwards_sol_matrix(rho_max_cut - rho_min_cut, backward_normalized_solutions)
        
        ## Find the Best Zeta From Zeta Array
        best_fit_index, best_fit_deltasq = best_fit_zeta(rho_match, forward_match_matrix, backward_match_matrix)
        
        ## Test the Best Fit Index:
        TEST_RESULT = zeta_interval_test(best_fit_index, zeta_array)
        
        if TEST_RESULT == False:
            raise Exception(r'The chosen interval for $\zeta$ was a bad guess. Choose a new interval')
            return 0
        
        ## Create a new min and max zeta:
        zeta_min = zeta_array[best_fit_index - 1]
        zeta_max = zeta_array[best_fit_index + 1]
        
        if SURPRESS_READOUT == False:
            ## Print Output
            print('Time Elapsed {} sec'.format(round(timeit.default_timer() - start_time), 4))
            print('The Best Guess ' + r'$\zeta$ for ' + 'interation {} is {} with a deltasq minimum of {}'.format(i, zeta_array[best_fit_index], best_fit_deltasq))
    
        best_fit_zeta_array[i] = zeta_array[best_fit_index]
        best_fit_deltasq_array[i] = best_fit_deltasq
        
    best_backwards_sol = backward_normalized_solutions[best_fit_index]
    best_forwards_sol = forward_normalized_solutions[best_fit_index]
    
    return best_fit_zeta_array, best_fit_deltasq_array, rhoB, rhoF, best_backwards_sol, best_forwards_sol 
