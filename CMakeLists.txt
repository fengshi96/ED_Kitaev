cmake_minimum_required(VERSION 3.16)
project(BitHamiltonian)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# Flags of sources, headers and libraries
set_directory_properties(PROPERTIES COMPILE_DEFINITIONS_DEBUG "_DEBUG")
#include_directories(/opt/OpenBLAS/include)
#link_directories(/opt/OpenBLAS/lib)
include_directories(
        /home/shifeng/Codes/0.Codes/intel/oneapi/mkl/latest/include
)

link_directories(
        /home/shifeng/Codes/0.Codes/intel/oneapi/mkl/latest/lib/intel64
        /home/shifeng/Codes/0.Codes/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin
)
set(MKL -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -liomp5 -lpthread)

set(SOURCE_FILES main.cpp src/Matrix.h src/Matrix.cpp src/Vector.h src/Vector.cpp src/LAPACK.h src/BLAS.h
        src/Parser.h src/QBasis.h src/QBasis.cpp src/Parser.cpp src/Hamiltonian.cpp
        src/Hamiltonian.h src/Lattice.cpp src/Lattice.h src/Memusage.h src/Memusage.cpp src/Lanczos.h)

# executable and linking
add_executable(BitHamiltonian ${SOURCE_FILES})
target_link_libraries(BitHamiltonian ${MKL})  #${MKL}