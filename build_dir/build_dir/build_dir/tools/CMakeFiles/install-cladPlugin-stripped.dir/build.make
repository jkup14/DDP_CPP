# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.27.6/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.27.6/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/clad

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir

# Utility rule file for install-cladPlugin-stripped.

# Include any custom commands dependencies for this target.
include tools/CMakeFiles/install-cladPlugin-stripped.dir/compiler_depend.make

# Include the progress variables for this target.
include tools/CMakeFiles/install-cladPlugin-stripped.dir/progress.make

tools/CMakeFiles/install-cladPlugin-stripped:
	cd /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/tools && /opt/homebrew/Cellar/cmake/3.27.6/bin/cmake -DCMAKE_INSTALL_COMPONENT="cladPlugin" -DCMAKE_INSTALL_DO_STRIP=1 -P /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/cmake_install.cmake

install-cladPlugin-stripped: tools/CMakeFiles/install-cladPlugin-stripped
install-cladPlugin-stripped: tools/CMakeFiles/install-cladPlugin-stripped.dir/build.make
.PHONY : install-cladPlugin-stripped

# Rule to build all files generated by this target.
tools/CMakeFiles/install-cladPlugin-stripped.dir/build: install-cladPlugin-stripped
.PHONY : tools/CMakeFiles/install-cladPlugin-stripped.dir/build

tools/CMakeFiles/install-cladPlugin-stripped.dir/clean:
	cd /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/tools && $(CMAKE_COMMAND) -P CMakeFiles/install-cladPlugin-stripped.dir/cmake_clean.cmake
.PHONY : tools/CMakeFiles/install-cladPlugin-stripped.dir/clean

tools/CMakeFiles/install-cladPlugin-stripped.dir/depend:
	cd /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/clad /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/clad/tools /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/tools /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/tools/CMakeFiles/install-cladPlugin-stripped.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : tools/CMakeFiles/install-cladPlugin-stripped.dir/depend

