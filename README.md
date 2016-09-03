# EPtran
EP transport code

This code can use OpenMP to parallize the main transport loop.

Two cases are included, 'ITERcaseAlpha' is the alpha particale in ITER 
and 'NBI_DIIIDcaseD' is the NBI (D partile) in DIII-D. 

The transport from micro-turbulence could use 
Angioni model (energy-dependent and energy-independent), 
DEP model (read input files), and Pueschel model.

The transport from AE is simulated by a critical gradient model, 
either critical density gradient or critical pressure gradient. 
The orbit drift effect can be included for this gradient.
C_R (F_cor) from TGLF simulations has been added to EPtran threshold. 
The new TGLF thresholds has also been tested.

The EPtran has been totally individual to the ALPHA code.

Most of the control parameters have been set as 'logical'.

The time-evolution uses the 2nd Runge-Kutta method.

The radial space and (E, lambda) velocity space are available.

The pitch angle (lambda) scattering is available.

The source (in E space) is still mono-energy and no time dependent.

The radial diffusion (D_rr) works well. 
The full 2x2 diffusivity matrix may have some numerical problems.
