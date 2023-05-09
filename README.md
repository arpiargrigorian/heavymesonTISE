# heavymesonTISE
Solves the Time Independent Schrodinger Equation with the Cornell Potential in three dimentions. This problem describes various heavy quark mesons, such as $c \bar{c}$, $b \bar{b}$, and $c \bar{b}$, where there are only two quarks and realtivistic effects can be ignored.  

## Parameters of the Cornell Pontential

The Cornell Potential is: 

$$ V(r) = \frac{-\kappa}{r} + \frac{r}{a^2} $$

The potential has a linear potential term and a Coulomb ($\propto 1/r$) term.

If the reduced mass of the system is $\mu$, we can write write a unitless form of the constant $\kappa$ as $\lambda$ with the relation:

$$ \lambda = \kappa (\mu a)^{2/3}$$

## Finding the Energy Spectrum as a Function of $\lambda$

### Inputs

To find the energy spectrum as a function of $\lambda$, *open the file inputs.py* and select your choice for the following parameters.

(Parentheses show the variable names of each parameter as they appear in the code).

_Constants of the Problem_

- The reduced mass of the two body system $\mu$ (`mu`) in units of GeV/c$^{2}$. The mass must be positive. 
- The constant for the strength of the linear potential $a$ (`a`) in units of (GeV/c$^{2}$)$^{-1}$. The constant $a$ must also be positive. 
- An array of numbers $\lambda$ (`coulomb_constant_array`) to test. $\lambda$ must be $geq$ 0.  $\lambda$ is unitless. The array of $\lambda$ can also contain one number. 
- A list of tuples that contains all the combinations of quantum numbers $n$, $l$ that will be in the spectrum. The first number in the tuple is the prinicpal quantum number $n$ and the second is the total angular momentum quantum number $l$. $n$ must be greater or equal to 1, $l$ is greater or equal to 0, and both $n$, $l$ must be integers. There is no need to order each tuple. 
  - i.e. `quantum_levels = [(1, 0), (2, 2), (2, 0)]`

_The Indepedent Variable_
- `rho`, `rho_cut_min`, `rho_cut_max`: 
  - `rho`: Is an array of ordered numbers representing the dimentionless distance $ \rho = a^{-1}(\mu a)^{1/3}
  - `rho_cut_min`: A positive integer. The value `rho[rho_cut_min]` is the smallest number that the *backwards* moving solution will solve to. 
  - `rho_cut_max`: A positive integer. The value `rho[rho_cut_max]` is the largest number that the *forwards* moving solution will solve to. 
Note: The Cornell + Centrigufal effective potential has singularities in both the backwards and forwards solutions. A good choice for `rho_cut_min` and `rho_cut_max` depends on the step size for `rho`. Be sure to test a few values of these parameters for a state $l > 0$ before solving over a larger parameter space. 

_The Eigenvalue Parameters_

- `zeta_init`: A dictionary of 2-tuples with the elements of the variable `quantum_levels` as keys. E.g.: `zeta_ini = {(1, 0): (.5, 2.5)}`. The tuple (a, b) represents the interval the code will search to find the best guess for the eigenvalue $\zeta$.
- `zeta_num_array`: A list or 1D array of positive integers. The code will iterate `zeta_num_array[0]` times first over the whole interval. Then, after finding the minimum $\zeta_{min}$, it will continue to iterate over a smaller interval ($\zeta_min - \delta \zeta, \zeta_min + \delta \zeta$) `zeta_num_array[1]` times. This process continues until the last element of `zeta_num_array`. A good guess for `zeta_num_array` is [40, 40, 40]. 

## Run Program

Once all inputs are specified, open the file _create_energy_spectrum.py_.

If you want to save outputs of the code as numpy files on your machine, change the variable SAVE_OUTPUT = True and specify the path you want to save the files in.

You can also choose to simply save the text output of the code if you run the file _create_energy_spectrum.py_ on the command line. This information only contains the best fit zeta value for each iteration and each combination of $n, l$ and $\lambda$. Just write:

`python create_energy_spectrum.py > simulation_out.txt`
