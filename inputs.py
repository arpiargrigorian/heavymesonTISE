## Import Statements
## Add new libraries as needed
import numpy as np

## Edit the Following Inputs as save before running create_energy_spectrum.py

mu = 1.85 ## Recommended Value is 1.85 GeV from Eichten+1978
a = 1.95 ## Recommended Value is 1.96 GeV^{-1} from Eichten+1978
coulomb_constant_array = np.arange(0, 1.5, num = 100)
levels = [(1, 0), (1,1), (1, 2), (2, 0), (2, 1), (2, 2)] ### n,l must be integers, n >= 1, l >= 0. 
