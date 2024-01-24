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

# Utility rule file for check-clad-enzyme.

# Include any custom commands dependencies for this target.
include test/CMakeFiles/check-clad-enzyme.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/check-clad-enzyme.dir/progress.make

test/CMakeFiles/check-clad-enzyme:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold --progress-dir=/Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Running lit suite /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/clad/test/Enzyme"
	cd /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/test && /Users/joshuakuperman/opt/anaconda3/envs/CPP/bin/lit -sv --param clad_site_config=/Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/test/lit.site.cfg /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/clad/test/Enzyme

check-clad-enzyme: test/CMakeFiles/check-clad-enzyme
check-clad-enzyme: test/CMakeFiles/check-clad-enzyme.dir/build.make
.PHONY : check-clad-enzyme

# Rule to build all files generated by this target.
test/CMakeFiles/check-clad-enzyme.dir/build: check-clad-enzyme
.PHONY : test/CMakeFiles/check-clad-enzyme.dir/build

test/CMakeFiles/check-clad-enzyme.dir/clean:
	cd /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/test && $(CMAKE_COMMAND) -P CMakeFiles/check-clad-enzyme.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/check-clad-enzyme.dir/clean

test/CMakeFiles/check-clad-enzyme.dir/depend:
	cd /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/clad /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/clad/test /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/test /Users/joshuakuperman/Desktop/Research/BaS_Cpp/build_dir/build_dir/build_dir/test/CMakeFiles/check-clad-enzyme.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : test/CMakeFiles/check-clad-enzyme.dir/depend

