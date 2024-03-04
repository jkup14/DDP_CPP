# Differential Dynamic Programming in C++
The popular trajectory optimization technique implemented in C++ and the Eigen Library for the linear algebra. Heavily typed and templated so that all matrices have knowns sizes. This allows Eigen to optimize during compilation. **Highlight: Sub 10 ms runtimes for Differential Drive with randomized start and goal states.**

**Next steps: additional dynamics functions, multi-shooting/multi-threading, mpc, automatic differentiation, control constraints via box-ddp**

For constrained ddp using barrier states and differentiable control, check out my python repo https://github.com/jkup14/Safe-PDP-DDP

Joshua Kuperman (jkup14@gmail.com, https://www.linkedin.com/in/joshuakuperman/)

Last Update 02/12/2024

## Dependencies (installed and updated as submodules automatically using CMake)
eigen 3.4.0
matplotplusplus 1.2.0


## Build instructions
~~~
cd /DDP_CPP

cmake -B build # make build directory

cmake --build build # make executables in /executables, run this every time you make a change
~~~

## Run Examples
~~~
./executables/Differential_Drive
./executables/Double_Integrator
./executables/Single_Integrator
./executables/Cart_Pole
~~~

## Code Structure

Example files **Double_Integrator.cpp, Single_Integrator.cpp, Differential_Drive.cpp, Cart_Pole.cpp** give a general idea on how the code should work.  

$~$

#### Headers in /include, .cpp's in /source

cost.hpp/cpp: cost interface and quadratic cost 

$~$

dynamics.hpp/cpp: dynamics interface

dynamics_derived: single and double integrator, differential drive, cart pole 

$~$

integrator.hpp/cpp: integrator interface, Euler integrator

$~$

ddp.hpp/cpp: DDP algorithm with linesearch and regularization implemented. 

$~$

plot_functions.hpp/cpp: Visualizer object with 2d visualization implemented 

