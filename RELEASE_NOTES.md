# MicroHH-DEPAC v1.0.0: Coupling of MicroHH Large Eddy Simulation with DEPAC Dry Deposition Module

## Citation
This release contains the source code used in:

> Rashidi, M., Krol, M. C., and van Zanten, M. C.: Coupling of MicroHH large eddy simulation with DEPAC dry deposition module, Geoscientific Model Development, in preparation, 2025.

## What this release contains
- Complete source code for MicroHH with DEPAC integration
- Fortran wrapper subroutine (`wrapper_depac.f90`) for C++/Fortran language coupling
- Modified `deposition.cxx` and `deposition.h` implementing the DEPAC-based NH₃ deposition scheme
- Modified `CMakeLists.txt` for mixed-language compilation
- Input configuration files (`.ini`) for all simulation cases (grassland, forest, heterogeneous)
- Post-processing Python scripts for figure generation

## Dependencies
- MicroHH base code: https://github.com/microhh/microhh
- DEPAC module: van Zanten et al. (2010)
- CMake ≥ 2.8.12
- gfortran
- NetCDF

## How to compile
See `README.md` for installation and compilation instructions.

## License
This code is distributed under the GNU General Public License v3.0.
