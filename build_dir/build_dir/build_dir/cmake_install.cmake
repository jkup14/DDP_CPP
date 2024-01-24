# Install script for directory: /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/clad

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/inst")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/clad/include/" FILES_MATCHING REGEX "/[^/]*\\.def$" REGEX "/[^/]*\\.h$" REGEX "/\\.svn$" EXCLUDE)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/clad/tools" TYPE DIRECTORY FILES "/Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/clad/tools/" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/include/" FILES_MATCHING REGEX "/cmakefiles$" EXCLUDE REGEX "/[^/]*\\.inc$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/lib/cmake_install.cmake")
  include("/Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/tools/cmake_install.cmake")
  include("/Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/test/cmake_install.cmake")
  include("/Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/demos/ErrorEstimation/CustomModel/cmake_install.cmake")
  include("/Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/demos/ErrorEstimation/PrintModel/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
