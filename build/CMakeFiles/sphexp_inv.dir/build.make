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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/david/fortran/sphexp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/david/fortran/sphexp/build

# Include any dependencies generated for this target.
include CMakeFiles/sphexp_inv.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/sphexp_inv.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/sphexp_inv.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/sphexp_inv.dir/flags.make

CMakeFiles/sphexp_inv.dir/src/sphexp_inv.f90.o: CMakeFiles/sphexp_inv.dir/flags.make
CMakeFiles/sphexp_inv.dir/src/sphexp_inv.f90.o: /home/david/fortran/sphexp/src/sphexp_inv.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/david/fortran/sphexp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/sphexp_inv.dir/src/sphexp_inv.f90.o"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/david/fortran/sphexp/src/sphexp_inv.f90 -o CMakeFiles/sphexp_inv.dir/src/sphexp_inv.f90.o

CMakeFiles/sphexp_inv.dir/src/sphexp_inv.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/sphexp_inv.dir/src/sphexp_inv.f90.i"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/david/fortran/sphexp/src/sphexp_inv.f90 > CMakeFiles/sphexp_inv.dir/src/sphexp_inv.f90.i

CMakeFiles/sphexp_inv.dir/src/sphexp_inv.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/sphexp_inv.dir/src/sphexp_inv.f90.s"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/david/fortran/sphexp/src/sphexp_inv.f90 -o CMakeFiles/sphexp_inv.dir/src/sphexp_inv.f90.s

# Object files for target sphexp_inv
sphexp_inv_OBJECTS = \
"CMakeFiles/sphexp_inv.dir/src/sphexp_inv.f90.o"

# External object files for target sphexp_inv
sphexp_inv_EXTERNAL_OBJECTS =

bin/sphexp_inv: CMakeFiles/sphexp_inv.dir/src/sphexp_inv.f90.o
bin/sphexp_inv: CMakeFiles/sphexp_inv.dir/build.make
bin/sphexp_inv: lib/libSphExpModules.a
bin/sphexp_inv: CMakeFiles/sphexp_inv.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/david/fortran/sphexp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable bin/sphexp_inv"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sphexp_inv.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/sphexp_inv.dir/build: bin/sphexp_inv
.PHONY : CMakeFiles/sphexp_inv.dir/build

CMakeFiles/sphexp_inv.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/sphexp_inv.dir/cmake_clean.cmake
.PHONY : CMakeFiles/sphexp_inv.dir/clean

CMakeFiles/sphexp_inv.dir/depend:
	cd /home/david/fortran/sphexp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/david/fortran/sphexp /home/david/fortran/sphexp /home/david/fortran/sphexp/build /home/david/fortran/sphexp/build /home/david/fortran/sphexp/build/CMakeFiles/sphexp_inv.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/sphexp_inv.dir/depend

