# LEGEND OF THE EXAMPLES



## Example 1

`ex1/`: Evaluation of the variational energy for a 1-particle 1-dimension harmonic oscillator using four different trial wave functions.



## Example 2

`ex2/`: Optimize a trial wave functions for -particle 1-dimension harmonic oscillator using the Conjugate Gradient method, setting the potential parameter w=1 and w=2.



## Example 3

`ex3/`: as `ex2`, but using a FFNN as wave function. It outputs some files that can be used to plot the wave function.
WARNING: Sometimes the conjugate gradient method get stuck and fails. Not always the optimisation method find a good minimum, depending on the starting variational parameters of the NN. (One could use always the same NN, but it can be instructive to see what happens with different starting points)
INFO: This example need the FFNN library.



## Example 4

`ex4/`: as `ex2`, but using the SR optimization method instead of the conjugate gradient



## Example 5

`ex5/`: as `ex2`, but using the Simulated Annealing method instead of the conjugate gradient



## Example 6

`ex6/`: Compute the variatonal energy of Gaussian vs. NN fitted to Gaussian, for 1D1P harmonic oscillator



## Example 7

`ex7/`: solve the 1 dimension 1 particle harmonic oscillator using a Neural Network Wave Function. The optimization is achieved with the Simulated Annealing algorithm.
INFO: a report folder with a collection of runs is stored inside this folder. For a quick overview on the results, just look at the file collage.jpg
