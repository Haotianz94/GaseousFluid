#GaseousFluid 
A Gaseous Fluid Simulation Demo

##Features
 * Based on the algorithm proposed in Jos Stam's paper "Real-time fluid dynamics for games"
 * Used Central Grid
 * No free surface presented
 * Used several methods for velocity field visualization
 * Used freeglut (alternative to OpenGL) for real-time rendering
 * Programed using C++

##Algorithm Framework

###For each iteration:

###Update velocity field based on the Navier-Strokes Equations
 * Apply external forces
 * Apply viscosity using Gauss-Seidel Relaxtion (diffusion)
 * Apply advection using "semi- Lagrangian" method
 * Apply the pressure projection using Gauss-Seidel Relaxtion
 * Set boundary conditions
 
###Update density field based on the Navier-Strokes Equations
 * Add density source
 * Apply viscosity using Gauss-Seidel Relaxtion (diffusion)
 * Apply advection using "semi- Lagrangian" method
 * Set boundary conditions
 
##References
[STAM, J., ¡°Real-time fluid dynamics for games,¡± in Proceedings of the Game Developer Conference, Mar. 2003.]
(http://www.intpowertechcorp.com/GDC03.pdf)
[Brian C., Leedom L., "Imaging Vector Fields using Line Integral Convolution," Proceedings of the 20th
annual conference on Computer graphics and interactive techniques. SIGGRAPH '93. Anaheim, California, 1993.]
(http://cs.brown.edu/courses/csci2370/2000/1999/cabral.pdf)

##Links
[Detailed description on my homepage](http://zhanghaotian1994.com/projects/GaseousFluid/)