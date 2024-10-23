cmake_minimum_required(VERSION 3.16) #last feature: target_precompile_header

# set the project name and version
project(NEGF VERSION 1.0)
set(CMAKE_VERBOSE_MAKEFILE ON)


option(NEGF_PROFILE "Profiling function timings" ON)
option(NEGF_MPI "Enabling MPI" ON)

if(NEGF_PROFILE)
  add_definitions(-DNEGF_PROFILE)
endif()
if(NEGF_MPI)
  add_definitions(-DNEGF_MPI)
  #set(MPI_CXX_COMPILER mpiicpc)
endif()


# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /lib/x86_64-linux-gnu/)
if (CMAKE_BUILD_TYPE="Debug")
    add_compile_options(-g -Wall)#-fconcepts -g)
endif()
add_compile_options(-fopenmp)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")



#git version
# check if git command exists
find_program(GIT_EXE NAMES git)

# generate version header
string(TIMESTAMP NEGF_TIMESTAMP "%Y-%m-%d %H:%M:%S")
message("PROJECT_SOURCE_DIR:     ${PROJECT_SOURCE_DIR}")
if(DEFINED GIT_EXE AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  execute_process(COMMAND git rev-parse HEAD
                  OUTPUT_VARIABLE NEGF_SHA
                  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
                  ERROR_QUIET
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND git describe --all
                  OUTPUT_VARIABLE NEGF_GIT_BRANCH
                  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
                  ERROR_QUIET
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "git hash ${NEGF_SHA}")
endif()

configure_file("${PROJECT_SOURCE_DIR}/src/core/githash.h.in" "${PROJECT_SOURCE_DIR}/src/core/githash.hpp" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/core/projectdir.h.in" "${PROJECT_SOURCE_DIR}/src/core/projectdir.hpp" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN.json.in" "${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN.json" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN_coulomb.json.in" "${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN_coulomb.json" @ONLY)

#compile sources as library INTERFACE (not a real library, just objects to be reused)
include_directories("${PROJECT_SOURCE_DIR}/src/")

set(_SOURCES src/initialize.cpp
             src/GlobalFunctions.cpp
             src/core/print_header.cpp
             src/ConvertUnits.cpp
             src/ostream.cpp
             src/StreamFile.cpp
             src/Geometry/Matrix.cpp 
             src/Geometry/Coordinate.cpp
             src/MeshGrid/MeshGrid.cpp
             src/Wannier/Wannier.cpp
             src/Wannier/ReadWannier.cpp
             src/Wannier/PrintWannier.cpp
             src/Laser/Laser.cpp
             src/fftPair/fftPair.cpp
             src/kGradient/kGradient.cpp
             src/RungeKutta/RungeKutta.cpp
             src/Simulation/Simulation.cpp
             src/core/profiler.cpp 
             src/core/print_timing.cpp 
             src/core/rt_graph.cpp
             src/core/mpi/Communicator.cpp
             src/Coulomb/Coulomb.cpp
             src/RytovaKeldysh/RytovaKeldysh.cpp
    )
	


#adding all sources to target
add_library(OBJECTS OBJECT ${_SOURCES} )
link_libraries(OBJECTS)

#adding OpenMP
find_package(OpenMP)
link_libraries(OpenMP::OpenMP_CXX)

#adding MPI
if(NEGF_MPI)
  find_package(MPI REQUIRED)
  link_libraries(MPI::MPI_CXX)
  message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
  if(MPI_COMPILE_FLAGS)
    set(COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    set(LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif()



#adding mkl (can it be replaced with lapack?)
set(MKL_ARCH intel64)
set(MKL_THREADING sequential)
set(MKL_INTERFACE lp64)
find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
add_compile_options($<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
include_directories($<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
link_libraries($<LINK_ONLY:MKL::MKL>)
message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")


#fftw3

#find_package(fftw3 REQUIRED)
link_libraries(fftw3)

if (NEGF_MPI)
  link_libraries(fftw3_mpi)
endif()

# add the executable
add_executable(EDUS src/EDUS.cpp)

#precompile headers to speed up compilation
#target_precompile_headers( NEGF 
#                              PUBLIC src/core/print_timing.hpp
#                              PUBLIC src/core/profiler.hpp
#                              PUBLIC src/core/rt_graph.hpp
#                              PUBLIC src/fftPair/fftPair.hpp
#                              PUBLIC src/Geometry/Coordinate.hpp
#                              PUBLIC src/Geometry/Matrix.hpp
#                              PUBLIC src/Geometry/Vector.hpp
#                              PUBLIC src/kGradient/kGradient.hpp
#                              PUBLIC src/Laser/Laser.hpp
#                              PUBLIC src/LinearAlgebra/gemm.hpp
#                              PUBLIC src/mdContainers/mdContainers.hpp
#                              PUBLIC src/MeshGrid/MeshGrid.hpp
#                              PUBLIC src/Model/Model.hpp
#                              PUBLIC src/Operator/BlockMatrix.hpp
#                              PUBLIC src/Operator/Operator.hpp
#                              PUBLIC src/RungeKutta/RungeKutta.hpp
#                              PUBLIC src/Simulation/Simulation.hpp
#                              PUBLIC src/Wannier/Wannier.hpp
#                              PUBLIC src/Constants.hpp
#                              PUBLIC src/ConvertUnits.hpp
#                              PUBLIC src/ostream.hpp
#                              PUBLIC src/StreamFile.hpp)




#testing...

#mdcontainers
add_executable(Containers_test ci-test/test_mdarray.cpp) 

#Matrix
add_executable(Matrix_test ci-test/test_matrix.cpp)

#Eigenvalues
add_executable(Eigen_test ci-test/test_eigen.cpp) 

#RK
add_executable(RK_test ci-test/test_RK.cpp)

#RK2
add_executable(RK2_test ci-test/test_RK2.cpp)

#fft-many
add_executable(fft-many_test ci-test/test_fft_many.cpp)

#fft+mpi
add_executable(fft_test ci-test/test_fft_mpi.cpp)

#Gradient
add_executable(Gradient_R_test ci-test/test_Gradient_R.cpp)

#Noneq propagation 
add_executable(NonEQ ci-test/test_NEQ_H0.cpp)

#Noneq propagation3 
add_executable(NonEQ3 ci-test/test_NEQ_Grad.cpp)

#bandstructure
add_executable(Bandstructure_test ci-test/test_bandstructure.cpp)

#pseudoinverse
add_executable(Pseudoinv_test ci-test/test_pseudoinverse.cpp)

#build for ctest
enable_testing()
add_test(NAME Containers_test COMMAND Containers_test)
add_test(NAME Matrix_test COMMAND Matrix_test)
add_test(NAME Eigen_test COMMAND Eigen_test)
add_test(NAME RK_test COMMAND RK_test)
add_test(NAME RK2_test COMMAND RK2_test)
add_test(NAME Gradient_R_test COMMAND Gradient_R_test)
add_test(NAME Bandstructure_test COMMAND 
					bash -c "${CMAKE_CURRENT_BINARY_DIR}/Bandstructure_test ; \
					python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_CURRENT_BINARY_DIR}/BANDSTRUCTURE.txt ${CMAKE_SOURCE_DIR}/ci-test/Reference/BANDSTRUCTURE.txt")
add_test(NAME fft_test COMMAND fft_test)
add_test(NAME NonEQ_H0 COMMAND NonEQ)
add_test(NAME NonEQ_Grad COMMAND NonEQ3)
add_test(NAME hBN COMMAND
             bash -c "${CMAKE_CURRENT_BINARY_DIR}/EDUS  ${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN.json; \
             python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py Population.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN/Population.txt
             python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py Velocity.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN/Velocity.txt; \
             python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py DM0.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN/DM0.txt")

add_test(NAME hBN_Coulomb COMMAND 
          bash -c "python3 ${PROJECT_SOURCE_DIR}/Postproces/RytovaKeldysh.py 40 40 1 ${PROJECT_SOURCE_DIR}/tb_models/hBN_gap7.25eV_a2.5A_tb.dat; \
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_CURRENT_BINARY_DIR}/RytovaKeldysh.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN_coulomb/RytovaKeldysh.txt; \
          ${CMAKE_CURRENT_BINARY_DIR}/EDUS  ${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN_coulomb.json; \
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_CURRENT_BINARY_DIR}/HF.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN_coulomb/HF.txt; \
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_CURRENT_BINARY_DIR}/Population.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN_coulomb/Population.txt; \
                  python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_CURRENT_BINARY_DIR}/Velocity.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN_coulomb/Velocity.txt")


