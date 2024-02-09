# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/joshuakuperman/Desktop/Research/DDP_CPP/build/_deps/matplotplusplus-src"
  "/Users/joshuakuperman/Desktop/Research/DDP_CPP/build/_deps/matplotplusplus-build"
  "/Users/joshuakuperman/Desktop/Research/DDP_CPP/build/_deps/matplotplusplus-subbuild/matplotplusplus-populate-prefix"
  "/Users/joshuakuperman/Desktop/Research/DDP_CPP/build/_deps/matplotplusplus-subbuild/matplotplusplus-populate-prefix/tmp"
  "/Users/joshuakuperman/Desktop/Research/DDP_CPP/build/_deps/matplotplusplus-subbuild/matplotplusplus-populate-prefix/src/matplotplusplus-populate-stamp"
  "/Users/joshuakuperman/Desktop/Research/DDP_CPP/build/_deps/matplotplusplus-subbuild/matplotplusplus-populate-prefix/src"
  "/Users/joshuakuperman/Desktop/Research/DDP_CPP/build/_deps/matplotplusplus-subbuild/matplotplusplus-populate-prefix/src/matplotplusplus-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/joshuakuperman/Desktop/Research/DDP_CPP/build/_deps/matplotplusplus-subbuild/matplotplusplus-populate-prefix/src/matplotplusplus-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/joshuakuperman/Desktop/Research/DDP_CPP/build/_deps/matplotplusplus-subbuild/matplotplusplus-populate-prefix/src/matplotplusplus-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
