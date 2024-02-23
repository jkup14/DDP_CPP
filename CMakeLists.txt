cmake_minimum_required(VERSION 3.12)
project(DDP_CPP CXX)
set(CMAKE_CXX_STANDARD 20)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ../executables)

add_subdirectory(matplotplusplus)


# FetchContent_Declare(matplotplusplus
#         GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus)

# FetchContent_GetProperties(matplotplusplus)
# if(NOT matplotplusplus_POPULATED)
#     FetchContent_Populate(matplotplusplus)
#     add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
# endif()

##### Single Integrator  #####
add_executable(Single_Integrator # Name of executable
    Single_Integrator.cpp # file with main in it, change this for new scripts
    # These shouldn't change. Hpps are built because of template linking errors
    include/cost.hpp
    include/ddp.hpp
    include/return_type_structs.hpp
    include/dynamics.hpp
    include/integrator.hpp
    source/plot_functions.cpp
    include/sign_func.hpp
    
)
target_include_directories(Single_Integrator #must be the same executable name
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    "/opt/homebrew/include/eigen3"
)
target_link_libraries(Single_Integrator PUBLIC
    matplot
)
###############################

##### Double Integrator  #####
add_executable(Double_Integrator
    Double_Integrator.cpp
    include/cost.hpp
    include/ddp.hpp
    include/return_type_structs.hpp
    include/dynamics.hpp
    include/integrator.hpp
    source/plot_functions.cpp
    include/sign_func.hpp
)
target_include_directories(Double_Integrator PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    "/opt/homebrew/include/eigen3"
)
target_link_libraries(Double_Integrator PUBLIC
    matplot
)
################################

##### Differential Drive  #####
add_executable(Differential_Drive
    Differential_Drive.cpp
    include/cost.hpp
    include/ddp.hpp
    include/return_type_structs.hpp
    include/dynamics.hpp
    include/integrator.hpp
    source/plot_functions.cpp
    include/sign_func.hpp
)
target_include_directories(Differential_Drive PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    "/opt/homebrew/include/eigen3"
)
target_link_libraries(Differential_Drive PUBLIC
    matplot
)
#########################


set(CMAKE_BUILD_TYPE Debug)  # Or Release