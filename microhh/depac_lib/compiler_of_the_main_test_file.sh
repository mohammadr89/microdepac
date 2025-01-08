#!/bin/bash

rm *.{mod,o}
# Compile Fortran modules
#   -cpp: Enable the C preprocessor, allowing for preprocessing directives.
#   -ffree-line-length-none: Allow for long lines in free form source code.
#   -c: This flag tells the compiler to compile and assemble, 
#       but do not link.
#  	  It compiles the source files into object (*.o) files
#
#
#   -fPIC: Stands for Position Independent Code. 
#          This is necessary when creating shared libraries.
#   .mod Files: contain information about Fortran modules, which are used to 
#  			  share variables, functions, and subroutines between different
# 			  parts of the program.

gfortran -c -fPIC -cpp -ffree-line-length-none go.f90
gfortran -c -fPIC -cpp -ffree-line-length-none le_drydepos_gas_depac.f90
gfortran -c -fPIC -cpp -ffree-line-length-none wrapper_depac.f90
gfortran -c -fPIC -cpp -ffree-line-length-none lai_wrapper.f90

# Create shared library
gfortran -shared -o libdepac.so go.o le_drydepos_gas_depac.o wrapper_depac.o lai_wrapper.o


# Compile C++ code: compile the main.cxx file into an object file (main.o),
#                   which can later be linked with other object files or 
# 					libraries to create an executable or a library.
# g++: This is the GNU C++ compiler. It's used to compile C++ source code.
# -I.: This flag adds the current directory (represented by the dot .) to the include search path.
# 	   The -I flag is used to specify additional directories to search for header files.
# 	   By including the current directory, the compiler will look for 
#	   any included headers in the same directory as the source file.
g++ -c -I. main.cxx lai_wrapper.h 


# Link everything
#   -lgfortran flag ensures that any necessary Fortran runtime support is included 
#   -L.: Adds the current directory (.) to the library search path
#   -ldepac: Links the depac library (searches for libdepac.a or libdepac.so)
g++ -o depac main.o lai_wrapper.o -L. -ldepac -lgfortran


# Set library path

# Adds current directory (.) to the list of directories
# where shared libraries are searched for during runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.

# Run the program
./depac




