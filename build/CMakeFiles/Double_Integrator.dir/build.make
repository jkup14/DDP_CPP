# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.28.2/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.28.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/joshuakuperman/Desktop/Research/DDP_CPP

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/joshuakuperman/Desktop/Research/DDP_CPP/build

# Include any dependencies generated for this target.
include CMakeFiles/Double_Integrator.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Double_Integrator.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Double_Integrator.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Double_Integrator.dir/flags.make

CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.o: CMakeFiles/Double_Integrator.dir/flags.make
CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.o: /Users/joshuakuperman/Desktop/Research/DDP_CPP/Double_Integrator.cpp
CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.o: CMakeFiles/Double_Integrator.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/joshuakuperman/Desktop/Research/DDP_CPP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.o -MF CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.o.d -o CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.o -c /Users/joshuakuperman/Desktop/Research/DDP_CPP/Double_Integrator.cpp

CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joshuakuperman/Desktop/Research/DDP_CPP/Double_Integrator.cpp > CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.i

CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joshuakuperman/Desktop/Research/DDP_CPP/Double_Integrator.cpp -o CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.s

CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.o: CMakeFiles/Double_Integrator.dir/flags.make
CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.o: /Users/joshuakuperman/Desktop/Research/DDP_CPP/source/plot_functions.cpp
CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.o: CMakeFiles/Double_Integrator.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/joshuakuperman/Desktop/Research/DDP_CPP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.o -MF CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.o.d -o CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.o -c /Users/joshuakuperman/Desktop/Research/DDP_CPP/source/plot_functions.cpp

CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joshuakuperman/Desktop/Research/DDP_CPP/source/plot_functions.cpp > CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.i

CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joshuakuperman/Desktop/Research/DDP_CPP/source/plot_functions.cpp -o CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.s

# Object files for target Double_Integrator
Double_Integrator_OBJECTS = \
"CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.o" \
"CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.o"

# External object files for target Double_Integrator
Double_Integrator_EXTERNAL_OBJECTS =

/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: CMakeFiles/Double_Integrator.dir/Double_Integrator.cpp.o
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: CMakeFiles/Double_Integrator.dir/source/plot_functions.cpp.o
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: CMakeFiles/Double_Integrator.dir/build.make
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: matplotplusplus/source/matplot/libmatplot.a
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /opt/homebrew/lib/libjpeg.dylib
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /opt/homebrew/lib/libtiff.dylib
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.2.sdk/usr/lib/libz.tbd
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /opt/homebrew/lib/libpng.dylib
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.2.sdk/usr/lib/libz.tbd
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /opt/homebrew/lib/libpng.dylib
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /opt/homebrew/lib/libfftw3.dylib
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /opt/homebrew/lib/libfftw3f.dylib
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /opt/homebrew/lib/libfftw3l.dylib
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /opt/homebrew/lib/libfftw3_threads.dylib
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /opt/homebrew/lib/libfftw3f_threads.dylib
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /opt/homebrew/lib/libfftw3l_threads.dylib
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /opt/homebrew/lib/libfftw3_omp.dylib
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /opt/homebrew/lib/libfftw3f_omp.dylib
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: /opt/homebrew/lib/libfftw3l_omp.dylib
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: matplotplusplus/source/3rd_party/libnodesoup.a
/Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator: CMakeFiles/Double_Integrator.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/joshuakuperman/Desktop/Research/DDP_CPP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable /Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Double_Integrator.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Double_Integrator.dir/build: /Users/joshuakuperman/Desktop/Research/DDP_CPP/executables/Double_Integrator
.PHONY : CMakeFiles/Double_Integrator.dir/build

CMakeFiles/Double_Integrator.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Double_Integrator.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Double_Integrator.dir/clean

CMakeFiles/Double_Integrator.dir/depend:
	cd /Users/joshuakuperman/Desktop/Research/DDP_CPP/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/joshuakuperman/Desktop/Research/DDP_CPP /Users/joshuakuperman/Desktop/Research/DDP_CPP /Users/joshuakuperman/Desktop/Research/DDP_CPP/build /Users/joshuakuperman/Desktop/Research/DDP_CPP/build /Users/joshuakuperman/Desktop/Research/DDP_CPP/build/CMakeFiles/Double_Integrator.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/Double_Integrator.dir/depend

