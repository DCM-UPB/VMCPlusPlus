# LEGEND OF THE EXAMPLES

Make sure the examples are compiled, by running `./build.sh` in the project root folder.
Execute an example by switching into one of the example folders and running `./run.sh`.
Note that the actual example executables reside inside the `build/examples/` folder under the project's root.

## Common

In the folder `common/` you will find a file named `ExampleFunctions.hpp` that contains some simple wave function
and Hamiltonian implementations that are shared between the examples.


## Basic Energy Evaluation

`ex_basic/`: Evaluation of the variational energy for a 1-particle 1-dimension harmonic oscillator using four different trial wave functions.


## Conjugate Gradient Optimization

`ex_cgopt/`: Optimize a trial wave functions for 1-particle 1-dimension harmonic oscillator using the Conjugate Gradient method.


## Stochastic Reconfiguration Optimization

`ex_sropt/`: as `ex_cgopt`, but using the SR gradient for Stochastic Gradient Descent optimization.


## Nelder-Mead Simplex Optimization

`ex_nmsopt/`: as `ex_cgopt`, but using the Nelder-Mead Simplex algorithm.


## Adam Optimization

`ex_adamopt/`: as `ex_cgopt`, but using the Adam algorithm.


## Basic Usage with MPI

`ex_mpi/`: Simple example of MPI-VMC usage. This example will only be compiled if USE_MPI=1 in config.sh (and MPI library found)!
