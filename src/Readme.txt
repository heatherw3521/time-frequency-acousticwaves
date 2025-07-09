This file includes ...

YOU NEED ON PATH needed:
* fmmlib2d: to run functions that call the fast multipole method
* sinctransform: to run the inverse Fourier transform (change matlab/sinc1d.m)
* finufft: called within sinctransform (nonuniform ffts)
* Chebfun: to run some boundary types
* cmocean: to plot with objectively better colors
* FFTW: to run files needed in other packages on path (fmm, sinctransform, finufft)
* Zetatrap2D: to run files related to quadrature and kernels formulations

experiment_driver is a test file that should run and produce a solution for different situations. In order to fully comprehend this, look at run_experiment_driver.

trapping_and_freqspace_images aims to help understand the trapping phenomenon in frequency space.
