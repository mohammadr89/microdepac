# 
#  MicroHH
#  Copyright (c) 2011-2020 Chiel van Heerwaarden
#  Copyright (c) 2011-2020 Thijs Heus
#  Copyright (c) 2014-2020 Bart van Stratum
# 
#  This file is part of MicroHH
# 
#  MicroHH is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  MicroHH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#

cmake_minimum_required (VERSION 2.8.12)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/config)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_CXX_STANDARD 17)

#----------------------
# DEPAC Integration START
#----------------------
# Define paths
set(DEPAC_LIB_DIR ${CMAKE_SOURCE_DIR}/depac_lib)
set(DEPAC_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

# Enable Fortran and set flags
enable_language(Fortran)

# Set Fortran compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffree-line-length-none -fPIC")

# Set C++ flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# Ensure shared libraries are built correctly
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(BUILD_SHARED_LIBS ON)

# Add the library path
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/depac_lib)
#----------------------
# DEPAC Integration END
#----------------------

# Set the precision of the build.
if(NOT USESP)
  set(USESP FALSE)
endif()

if(USESP)
  add_compile_definitions(FLOAT_SINGLE)
  add_compile_definitions(RTE_RRTMGP_SINGLE_PRECISION)
  message(STATUS "Precision: Single (32-bits floats)")
else()
  message(STATUS "Precision: Double (64-bits floats)")
endif()

# Check whether USEMPI and USECUDA are set, if not, set to FALSE.
if(NOT USEMPI)
  set(USEMPI FALSE)
endif()
if(NOT USECUDA)
  set(USECUDA FALSE)
endif()

# Crash on using CUDA and MPI together, not implemented yet.
if(USEMPI AND USECUDA)
  message(FATAL_ERROR "MPI support for CUDA runs is not supported yet")
endif()

# Load system specific settings if not set, force default.cmake.
if(NOT SYST)
  set(SYST default)
endif()
include(${SYST} OPTIONAL RESULT_VARIABLE SYSTINC)

# Trigger fatal error if illegal module is loaded.
if(${SYSTINC} STREQUAL "NOTFOUND")
  message(FATAL_ERROR "Config file config/" ${SYST} ".cmake does not exist.")
endif()

# Set the default build type to RELEASE.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: None Debug Release." FORCE)
else()
  string(TOUPPER ${CMAKE_BUILD_TYPE} TEMP)
  set(CMAKE_BUILD_TYPE ${TEMP} CACHE STRING
    "Choose the type of build, options are: None Debug Release." FORCE)
endif()

# Start the project only after the system specific settings are loaded.
if(USECUDA)
  project(microhhc C CXX CUDA Fortran)
else()
  project(microhhc C CXX Fortran)
endif()



# Display status messages for MPI set precompiler flags.
if(USEMPI)
  message(STATUS "MPI: Enabled.")
  # add the necessary precompiler flag
  add_definitions("-DUSEMPI")
else()
  message(STATUS "MPI: Disabled.")
endif()

# Load the CUDA module in case CUDA is enabled and display status message.
if(USECUDA)
  message(STATUS "CUDA: Enabled.")
  add_definitions("-DUSECUDA")
else()
  message(STATUS "CUDA: Disabled.")
endif()

# Only set the compiler flags when the cache is created
# to enable editing of the flags in the CMakeCache.txt file.
if(NOT HASCACHE)
  set(CMAKE_C_FLAGS ${USER_C_FLAGS} CACHE STRING
    "Flags used by the C-compiler during all build types." FORCE)
  set(CMAKE_C_FLAGS_DEBUG ${USER_C_FLAGS_DEBUG} CACHE STRING
    "Flags used by the C-compiler during debug builds." FORCE)
  set(CMAKE_C_FLAGS_RELEASE ${USER_C_FLAGS_RELEASE} CACHE STRING
    "Flags used by the C-compiler during release builds." FORCE)
  set(CMAKE_CXX_FLAGS ${USER_CXX_FLAGS} CACHE STRING
    "Flags used by the CXX-compiler during all build types." FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG ${USER_CXX_FLAGS_DEBUG} CACHE STRING
    "Flags used by the CXX-compiler during debug builds." FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE ${USER_CXX_FLAGS_RELEASE} CACHE STRING
    "Flags used by the CXX-compiler during release builds." FORCE)
# Commented the original line
## set(CMAKE_Fortran_FLAGS ${USER_FC_FLAGS} CACHE STRING
##    "Flags used by the Fortran-compiler during all build types." FORCE)

# Added the preprocessor and line length options
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffree-line-length-none" CACHE STRING
    "Flags used by the Fortran-compiler during all build types." FORCE)

  set(CMAKE_Fortran_FLAGS_DEBUG ${USER_FC_FLAGS_DEBUG} CACHE STRING
    "Flags used by the Fortran-compiler during debug builds." FORCE)
  set(CMAKE_Fortran_FLAGS_RELEASE ${USER_FC_FLAGS_RELEASE} CACHE STRING
    "Flags used by the Fortran-compiler during release builds." FORCE)

  if(USECUDA)
    set(CMAKE_CUDA_FLAGS ${USER_CUDA_NVCC_FLAGS} CACHE STRING
      "Flags used by NVCC during all build types." FORCE)
    set(CMAKE_CUDA_FLAGS_DEBUG ${USER_CUDA_NVCC_FLAGS_DEBUG} CACHE STRING
      "Flags used by NVCC during debug builds." FORCE)
    set(CMAKE_CUDA_FLAGS_RELEASE ${USER_CUDA_NVCC_FLAGS_RELEASE} CACHE STRING
      "Flags used by NVCC during release builds." FORCE)
  endif()

  message(STATUS "Build Type: " ${CMAKE_BUILD_TYPE})
  set(HASCACHE TRUE CACHE BOOL "CMakeCache.txt created." FORCE)
  # Make sure that ccmake only contains build type.
  mark_as_advanced(HASCACHE)
  mark_as_advanced(CMAKE_INSTALL_PREFIX)
endif()

# Print the C++ and CUDA compiler flags to the screen.
if(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  message(STATUS "Compiler flags: " ${CMAKE_CXX_FLAGS} " " ${CMAKE_CXX_FLAGS_RELEASE})
else()
  message(STATUS "Compiler flags: " ${CMAKE_CXX_FLAGS} " " ${CMAKE_CXX_FLAGS_DEBUG})
endif()

if(USECUDA)
  # Add CUDA debug definition, to enforce CUDA API and kernel checks.
  if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
   # add_definitions("-DCUDACHECKS")
    message(STATUS "NVCC flags: " ${CMAKE_CUDA_FLAGS} " " ${CMAKE_CUDA_FLAGS_DEBUG})
  else()
    message(STATUS "NVCC flags: " ${CMAKE_CUDA_FLAGS} " " ${CMAKE_CUDA_FLAGS_RELEASE})
  endif()
endif()

add_subdirectory(rte-rrtmgp-cpp/src_kernels)
if(USECUDA)
  add_subdirectory(rte-rrtmgp-cpp/src_kernels_cuda)
endif()
add_subdirectory(rte-rrtmgp-cpp/src)
if(USECUDA)
  add_subdirectory(rte-rrtmgp-cpp/src_cuda)
endif()
add_subdirectory(src)
add_subdirectory(main)
add_subdirectory(kpp)


