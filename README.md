# EPtran
EP transport code

Two basic cases are introduced, 
'ITERcaseAlpha' is the alpha particales in ITER with D-T background plasma
and 'NBI_DIIIDcaseD' is the NBI (D particles) in DIII-D with pure D background plasma. 

The transport from micro-turbulence could be simulated by 
Angioni model (energy-dependent or energy-independent, default), 
DEP model (read from input files, with DEP simulations), 
and Pueschel model (with/without electromagnetic part).

The transport from Alfven eigenmode (AE) is simulated by a critical gradient model, 
either critical density gradient (default) or critical pressure gradient. 
The critical gradient model can be energy-indenpent (default) or energy-denpendent (G(E)).
The orbit drift effect can be included for this gradient.
C_R (F_cor) from TGLF simulations has been added to EPtran threshold (only for critical density gradint). 
The new TGLF thresholds has also been tested.

The EPtran has been totally individual to the ALPHA code.

Most of the control parameters have been set as 'logical' instead of 'int' in ALPHA.

This code can use OpenMP to parallize the main transport loop.

The time-evolution of transport equations uses the 2nd Runge-Kutta method.

The radial space and (E, lambda) velocity space are available.

The pitch angle (lambda) scattering is available.

The source (in E space) is still mono-energy and time independent.

The radial diffusion (D_rr) works well. 
The full 2x2 diffusivity matrix may have some numerical problems.
