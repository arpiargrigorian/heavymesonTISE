'''
Creates the Energy Spectrum as a function of the energy levels and lambda. 

You can either choose:
Save the output of the script as a text file. This will show the zeta values for each level. 
Choose SAVE_OUTPUT = True. This will write the solutions, their corresponding rho arrays, and the zeta_values/mind distances as a function of lambda to the path on your machine you specify with save_path_folder. 

'''

from inputs import *
from find_best_zeta_library import *

## Optional: Save path Folder.  
save_path_folder = '' #(add the / at the end of the path)

SAVE_OUTPUT = False # Make this True if you want the simulation results saved. 

for quantum_level in levels:
    
    zeta_function = np.zeros(coulomb_constant_array.shape)
    distance_function = np.zeros(coulomb_constant_array.shape)
    
    for i, coulomb_constant in enumerate(coulomb_constant_array):
        
        el = quantum_level[1]
        quantum_n = quantum_level[0]

        rho_cut_interval = (rho_cut_min, rho_cut_max)

        zeta_init_interval = zeta_init[quantum_level]
        
        print('The Simulation Results for n = {}, l = {} and lambda = {} are:'.format(quantum_n, el, coulomb_constant))
        
        simulation_results = iter_find_zeta(el, quantum_n, coulomb_constant, rho, rho_cut_interval, zeta_init_interval, zeta_num_array)
        
        if SAVE_OUTPUT == True:
            
            np.save(save_path_folder + 'hB_n{}el{}_lambda{}.npy'.format(quantum_n, el,coulomb_constant), simulation_results[2])
            np.save(save_path_folder + 'hF_n{}el{}_lambda{}.npy'.format(quantum_n, el,coulomb_constant), simulation_results[3])
            
            np.save(save_path_folder + 'solB_n{}el{}_lambda{}.npy'.format(quantum_n, el,coulomb_constant), simulation_results[4])
            np.save(save_path_folder + 'solF_n{}el{}_lambda{}.npy'.format(quantum_n, el,coulomb_constant), simulation_results[5])
        
        zeta_function[i] = simulation_results[0][-1]
        distance_function[i] = simulation_results[1][-1]   
    
    if SAVE_OUTPUT == True:
        
        np.save(save_path_folder + 'zeta_function_n{}el{}.npy'.format(quantum_n, el), zeta_function)
        np.save(save_path_folder + 'distance_function_n{}el{}.npy'.format(quantum_n, el), distance_function)
