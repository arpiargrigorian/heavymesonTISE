## Import Statements
## Add new libraries as needed
import numpy as np

## Edit the Following Inputs as save before running create_energy_spectrum.py
## Run to check if all chosen parameters are valid. 

# Parameters of the Problem
mu = 1.85 ## Recommended Value is 1.85 GeV for charmonium (c cbar) from Eichten+1978
a = 1.95 ## Recommended Value is 1.96 GeV^{-1} for charmonium (c cbar) from Eichten+1978
coulomb_constant_array = np.linspace(0, 1.5, num = 3) # Numbers must be >= 0. For the charmonium problem, it is good to go up to 1.5. For the bottomonium problem, try a max value of 3.5. 
levels = [(1, 0), (1,1), (1, 2), (2, 0), (2, 1), (2, 2)] ### n,l must be integers, n >= 1, l >= 0.

# Indepdent Variable Parameters
rho = np.linspace(0, 30, num = 1000)
rho_cut_min = 5 #Just use a positive index here to avoid zero.
rho_cut_max = 250 ## You can use positive or negative values here. 

# NOTE: The rho values above avoid the singularities for every energy level EXCEPT the n = 1, l = 0 level. 
# For that level, you can run it seperatley with rho_cut_min = 5 and rho_cut_max = 250

# Eigenvalue Paramters
zeta_init = {(1, 0):(.5, 2.5), (1,1):(2.5, 3.5), (1, 2):(3.5, 4.5), (2, 0):(3, 4.5), (2, 1):(4, 5), (2, 2):(5, 6)} #The code will tell you if a guess is bad. Consult the measured spectrum for good guesses. 
zeta_num_array = [40, 40, 40]

# =========== TESTING PARAMTERS =================

assert mu > 0
assert a > 0
assert all(number >= 0 for number in coulomb_constant_array)

for quantum_numbers in levels:
    assert type(quantum_numbers[0]) == int and quantum_numbers[0] >= 1
    assert type(quantum_numbers[1]) == int and quantum_numbers[1] >= 0
  
assert all(a <= b for a, b in zip(rho, rho[1:]))
assert type(rho_cut_min) == int and rho_cut_min > 0
assert type(rho_cut_max) == int

assert type(zeta_init) == dict

for quantum_numbers in levels:
    assert quantum_numbers in zeta_init.keys() 


