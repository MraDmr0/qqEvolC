# qqEvolC (Ver. 1.0)
Package to simulate the time evolution of a D-level quantum system under the influence of a scalar potential.

## Description
This package implements a C++ version of the Runge-Kutta fourth order algorithm to integrate the time-dependent Schrodinger equation. 
In particular, this is applied to the determination of the time evolution of the occupation probability $|\Psi_j(t)|^2$ of a generic D-levels quantum system under the influence of an external oscillating scalar potential. The physical system is specified trough the position of its energy levels (Larmor frequencies $\omega_L$). The scalar potential is characterized by its oscillating frequency $\omega$ and a time-dependent envelope function $F(t)$ that allows to shape the potential. The coupling between potential and system is instead given by the Rabi frequencies $\omega_{i,j}$, where the indices $i,j$ label the energy levels. The overall potential that enters the Schrodinger equation reads:
```math
V(t)_{i,j} = -i F(t) cos( \omega t ) \omega _{i,j} \exp\left(-i( \omega ^ {j}_{L}-\omega^{i}_{L})t\right).
```
The current version supports the following envelope functions (for further details see "Input description"):
* "Off": vanishing potential
* "Constant": the envelope function is constant during the whole evolution
* "Single impulse": the envelope has the shape of one rectangular impulse
* "Double impulse": the envelope has the shape of two rectangular impulse
* "Single gaussian": the envelope has the shape of one gaussian impulses
* "Double gaussian": the envelope has the shape of two gaussian impulses

## Use
The source code can be compiled trough the existing Makefile by running the command:

$ make.

This will produce the executable "main" that than can be used as:

$ ./main.exe input.json > output.

! Note: currently, the Makefile only supports the GNU (g++) compiler. Different compilers can be easily integrated in the existing Makefile. 

## Input description
The input parameters can be specified in a JSON file. You need to specify a set of mandatory data needed to charaterize the basic system that you want to simulate, while the need for optional parameters depend on the choiche of the envelope function.

The input parameters can be specified in a JSON file. You need to specify a set of mandatory data needed to charaterize the basic system that you want to simulate, while the need for optional parameters depend on the choiche of the envelope function.

### Mandatory parameters
* prefix                = string used to save output
* D                     = number of energy levels
* ti                    = initial time ($\mu s$)
* tf                    = final time ($\mu s$)
* N                     = number of time steps   
* S                     = number of interations after which the result is saved
* [pis_0, ...]          = initial state
* [wl_0, ...]           = Larmor frequencies expressed in energy units (meV)
* [w_00, ...]           = Rabi frequencies of the system (Hz)
* qb_mode               = currently only supported the "off" mode
* env_mode              = specifies the envelope function. Supported modes:
    * "off"               : no potential
    * "const"             : constant potential in [ti , tf] with frequency w1 and amplitude F1
    * "impulse"           : one sqare impulse of frequency "w1" in [t0 , t1] of amplitude F1
    * "gauss"             : one gaussian impulse centered in "t0" of frequency "w1", spreading "sigma1" and amplitude "F1"
    * "double_impulse"    : two square impulses in [t0 , t1] (with frequency "w1" and amplitude "F1") and another in [t00 , t11] (with frewuency "w2" and amplitude "F2")
    * "double_gauss"      : two gaussian impulses centered in "t0" (with frequency "w1", spreading "sigma1" and amplitude "F1") and another in "t1" (with frequency "w2", spreading "sigma2" and amplitude "F2")
 
#### Optional parameters
* w1               = frequency of the first (or only) potential impulse (Hz)
* w2               = frequency of the second potential impulse (Hz)
* F1               = scaling fact of the first (or only) potential impulse
* F2               = scaling factor of the second potential impulse
* t0               = begin of first impulse if square, or center of the first gaussian impulse ($\mu s$)
* t1               = end of first impulse if square, or center of the second gaussian impulse ($\mu s$)
* t00              = begin of the second square impulse ($\mu s$)
* t11              = end of the second square impulse ($\mu s$)
* sigma1           = spreading of the first gaussian impulse
* sigma2           = spreading of the second gaussian impulse

## Output description
The executable outputs a file "prefix.txt" which contains, at each time step saved, the following data:
* time ($\mu s$)
* value of envelope function
* occupation probability of each energy level

## Utilities
The package includes a Python script "Plot.py" that allows to easily plot the occupation probability and the envelope shape as a function of time:

$python Plot.py prefix.

## Next developements
The current Ver. 1.0 implements the basic functions of the package. However, the code and the repository can be further developed from many points of view:
* Overall the code should be better organized and commentend to become more user and developer friendly
* A reorganization of the code should also aim to an enhancement of the performances, even with multi-thread parallelization if worth
* The code can be extended by implementing more envelope functions or different potentials
* The repository lacks practical examples which need to be included
* The python script is simple and only works up to four levels, this needs to be generalized
* The code was initially thought to implement also "qb_mode = on" to accept input parameters in a different form for the special case of a two-levels system

## Other informations
This code follows the Python version "Evoluzione.py" (https://github.com/MraDmr0/Evoluzione), which however lacks input and error handling. Moreover, the Python implementation is sevearly constrained by the computational performances offered by the language. Another Python version which leverages the Numba package for performance enhancements is under development.

Mario Di Mare, 01/10/2025
