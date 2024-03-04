cmake_minimum_required(VERSION 3.12)
project(DDP_CPP CXX)
set(CMAKE_CXX_STANDARD 20)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ../executables)

# Update the submodules here
include(cmake/UpdateSubmodules.cmake)

add_subdirectory(matplotplusplus)

# Add additional executables to this list
list(APPEND executables Single_Integrator Double_Integrator Differential_Drive Single_Integrator_test Cart_Pole)
# list(APPEND executables Single_Integrator_test)

foreach(exec ${executables})
    add_executable(${exec} # Name of executable
        "${exec}.cpp" 
        # These shouldn't change. Hpps are built because of template linking errors
        include/cost.hpp
        include/ddp.hpp
        include/return_type_structs.hpp
        include/integrator.hpp
        include/dynamics.hpp
        include/dynamics_derived/SingleIntegrator.hpp
        include/dynamics_derived/DoubleIntegrator.hpp
        include/dynamics_derived/DifferentialDrive.hpp
        include/dynamics_derived/CartPole.hpp
        source/timer.cpp
        source/plot_functions.cpp
        include/helper_funcs.hpp
    )
    target_include_directories(${exec} #must be the same executable name
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        eigen
        matplotplusplus/source
    )
    target_link_libraries(${exec} PUBLIC
        matplot
    )
endforeach()

set(CMAKE_BUILD_TYPE Release)  # Debug Or Release
