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

- The reduced mass of the two body system $\mu$ (`mu`) in units of GeV/c$^{2}$. The mass must be positive. 
- The constant for the stregth of the linear potential $a$ (`a`) in units of (GeV/c$^{2}$)$^{-1}$. 
- An array of $\lambda$ (`coulomb_constant_array`) to test. $\lambda$ is unitless. 
- A list of tuples that contains all the combinations of quantum numbers $n$, $l$ that will be in the spectrum. The first number in the tuple is the prinicpal quantum number $n$ and the second is the total angular momentum quantum number $l$. $n$ must be greater or equal to 1, $l$ is greater or equal to 0, and both $n$, $l$ must be integers. There is no need to order each tuple. 
  - i.e. `quantum_levels = [(1, 0), (2, 2), (2, 0)]`

(Later, I need to put in how other things are organized.) the file to run is create_energy_spectrum.py

For a two-body heavy-quark meson system, the linear term $\frac{r}{a^2}$ is t
