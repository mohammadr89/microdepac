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
# First create DEPAC library
set(DEPAC_SOURCES
    ${CMAKE_SOURCE_DIR}/depac_lib/go.f90
    ${CMAKE_SOURCE_DIR}/depac_lib/le_drydepos_gas_depac.f90
    ${CMAKE_SOURCE_DIR}/depac_lib/wrapper_depac.f90
)

# Create DEPAC library
add_library(depac STATIC ${DEPAC_SOURCES})

# Set DEPAC properties
set_target_properties(depac PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/depac_lib
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/depac_lib
    COMPILE_FLAGS "-cpp -ffree-line-length-none -fPIC"
)

# Original MicroHH configuration
FILE(GLOB sourcefiles "../src/*.cxx" "../src/*.cu")

include_directories(
    "../include" 
    "../rte-rrtmgp-cpp/include" 
    ${CMAKE_BINARY_DIR}/modules  # For Fortran modules
    SYSTEM ${INCLUDE_DIRS}
)

add_library(microhhc STATIC ${sourcefiles})

# Add dependency on depac
add_dependencies(microhhc depac)

target_include_directories(microhhc PRIVATE 
    "../include" 
    "../rte-rrtmgp-cpp/include" 
    "../rte-rrtmgp-cpp/include_rt" 
    "../rte-rrtmgp-cpp/include_kernels_cuda" 
    "../rte-rrtmgp-cpp/include_rt_kernels" 
    ${CMAKE_BINARY_DIR}/modules
    ${INCLUDE_DIRS}
)

if(USECUDA)
    if(USEKERNELLAUNCHER)
        target_link_libraries(microhhc PRIVATE CUDA::curand CUDA::cufft ${LIBS} kernel_launcher depac gfortran)
    else()
        target_link_libraries(microhhc PRIVATE CUDA::curand CUDA::cufft ${LIBS} depac gfortran)
    endif()
else()
    target_link_libraries(microhhc PRIVATE ${LIBS} depac gfortran)
endif()


