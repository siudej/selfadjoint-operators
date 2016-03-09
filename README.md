# Spectra of self-adjoint operators with four eigenvalues.

This Mathematica script is an optimized version of the code used by Marcin Bownik and John Jasper in their paper:
##### [Spectra of Frame Operators with Prescribed Frame Norms](http://pages.uoregon.edu/mbownik/papers/53.pdf)

The original code (`old.nb`) involves depth-first approach through various nested mathematical definitions and a nonlinear optimization step. 

The optimal code precomputes values of some functions, eliminates optimization step, and tabularizes all repeating computations.

As a result:
* constraint on multiplicity of eigenvalues is now eliminated
* all computations can be performed using rational numbers
* execution time decreased more than 1000 times. Finding images for {b, 0.752, 0.766, 0.002} (8 values):
  - old: 8 values computed in parallel, 1 threads each - 1500 s
  - new: 8 threads per value, time to finish all in serial - 0.2 s  


