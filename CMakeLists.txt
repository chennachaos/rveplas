cmake_minimum_required (VERSION 3.0)

project (hypas)

enable_language(Fortran)

##################################
# compiler options
##################################


set(CMAKE_Fortran_COMPILER "gfortran")

#debug
set(FFLAGS "-g -pg -O2 -fopenmp -cpp")


set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FLAGS}")


##################################
# include directories
##################################
#include_directories(/usr/include/petsc)
include_directories(
)


##################################
# lib directories
##################################
#link_directories(/usr/lib/gcc/x86_64-linux-gnu /usr/lib /opt/petsc-3.6.4/arch-linux2-c-debug/lib)
link_directories(
/usr/lib/gcc/x86_64-linux-gnu/7
/usr/lib/x86_64-linux-gnu 
/home/chenna/Documents/otherCodes/hyplas/pardisolib
/usr/lib/x86_64-linux-gnu/blas
/usr/lib/x86_64-linux-gnu/lapack
)



#src/elemutilitiesquadrature.F 
file(GLOB_RECURSE sources1 src/*/*.f)
file(GLOB_RECURSE sources2 src/*/*/*.f)

##################################

add_executable(rveplas src/hyplas.f ${sources1}  ${sources2})

target_link_libraries(rveplas gfortran gomp -fopenmp pardiso600-GNU720-X86-64 blas lapack ) # mkl_blacs_lp64 mkl_blas95_lp64 mkl_lapack95_lp64 mkl_scalapack_lp64 mkl_scalapack_ilp64 mkl_lapack95_lp64 mkl_lapack95_ilp64 mkl_intel_thread mkl_intel_ilp64 mkl_gnu_thread mkl_gf_lp64 mkl_gf_ilp64 mkl_blas95_ilp64 mkl_cdft_core mkl_blacs_openmpi_lp64 mkl_blacs_sgimpt_ilp64 mkl_blacs_sgimpt_lp64 mkl_blacs_intelmpi_lp64 mkl_blacs_openmpi_ilp64 mkl_blacs_ilp64 mkl_blacs_intelmpi_ilp64 )

install(TARGETS rveplas RUNTIME DESTINATION /home/chenna/Documents/otherCodes/hyplas/bin)
