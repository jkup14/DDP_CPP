# Differential Dynamic Programming in C++
The popular trajectory optimization technique implemented in C++ and the Eigen Library for the linear algebra. Heavily typed and templated so that all matrices have knowns sizes. This allows Eigen to optimize during compilation. 

**Next steps: visualization, additional dynamics functions, regularization, automatic differentiation, allocate on the heap, control constraints via box-ddp**

For constrained ddp using barrier states and differentiable control, check out my python repo https://github.com/jkup14/Safe-PDP-DDP

Joshua Kuperman (jkup14@gmail.com, https://www.linkedin.com/in/joshuakuperman/)

Last Update 02/09/2024

## Build instructions
~~~
cmake -B build # make build directory

cmake --build build # make executables in /executables, run this every time you make a change
~~~

## Code Structure

Example files **Double_Integrator.cpp, Single_Integrator.cpp, Differential_Drive.cpp** gives a general idea on how the code should work.  

$~$

#### Headers in /include, .cpp's in /source

cost.hpp: cost interface

cost.cpp: quadratic cost function  

$~$

dynamics.hpp: dynamics interface

dynamics.cpp: single and double integrator, differential drive  

$~$

integrator.hpp: integrator interface

integrator.cpp: Euler integrator for discretized dynamics  

$~$

ddp.hpp/cpp: DDP algorothm with linesearch implemented. Regularization will be added soon. 

