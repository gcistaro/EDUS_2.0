cmake_minimum_required(VERSION 3.16) #last feature: target_precompile_header

# set the project name and version
project(EDUS VERSION 1.0)
set(CMAKE_VERBOSE_MAKEFILE ON)


option(EDUS_PROFILE     "Profiling function timings"                            ON)
option(EDUS_MPI         "Enabling MPI"                                          ON)
option(EDUS_BATCHGEMM   "Enabling Batch Gemm for k matrix multiplication"      OFF)
option(EDUS_FFTWTHREADS "Enabling threads in fftw"                             OFF)                
option(EDUS_MKL_THREAD  "Enabling threads in mkl"                              OFF)                  
option(EDUS_HDF5        "Link hdf5 library for fast writing of big matrices"    ON)

#set(HDF5_PREFER_PARALLEL  true)
if(EDUS_PROFILE)
  add_definitions(-DEDUS_PROFILE)
endif()
if(EDUS_MPI)
  add_definitions(-DEDUS_MPI)
endif()
if(EDUS_BATCHGEMM)
  add_definitions(-DEDUS_BATCHGEMM)
endif()
if(EDUS_FFTWTHREADS)
  add_definitions(-DEDUS_FFTWTHREADS)
endif()
if(EDUS_MKL_THREAD)
  add_definitions(-DEDUS_MKL_THREAD)
endif()
if(EDUS_HDF5)
	add_definitions(-DEDUS_HDF5)
endif()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /lib/x86_64-linux-gnu/)

if (CMAKE_BUILD_TYPE="Debug")
    add_compile_options(-O0 -g -Wall)
endif()

if (CMAKE_BUILD_TYPE="Release")
    add_compile_options(-O3)
endif()



list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")



#git version
# check if git command exists
find_program(GIT_EXE NAMES git)

# generate version header
string(TIMESTAMP EDUS_TIMESTAMP "%Y-%m-%d %H:%M:%S")
message("PROJECT_SOURCE_DIR:     ${PROJECT_SOURCE_DIR}")
if(DEFINED GIT_EXE AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  execute_process(COMMAND git rev-parse HEAD
                  OUTPUT_VARIABLE EDUS_SHA
                  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
                  ERROR_QUIET
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND git describe --all
                  OUTPUT_VARIABLE EDUS_GIT_BRANCH
                  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
                  ERROR_QUIET
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "git hash ${EDUS_SHA}")
endif()

configure_file("${PROJECT_SOURCE_DIR}/src/core/githash.h.in" "${PROJECT_SOURCE_DIR}/src/core/githash.hpp" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/core/projectdir.h.in" "${PROJECT_SOURCE_DIR}/src/core/projectdir.hpp" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN.json.in" "${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN.json" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN_coulomb.json.in" "${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN_coulomb.json" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/ci-test/inputs/MoS2.json.in" "${PROJECT_SOURCE_DIR}/ci-test/inputs/MoS2.json" @ONLY)

include_directories("${PROJECT_SOURCE_DIR}/src/")
########add_subdirectory(src/)


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
             src/DESolver/DESolver.cpp
             src/Simulation/Simulation.cpp
             src/core/profiler.cpp
             src/core/print_timing.cpp
             src/core/cmd_args/cmd_args.cpp
             src/core/rt_graph.cpp
             src/core/mpi/Communicator.cpp
             src/Coulomb/Coulomb.cpp
             src/GreenFunction/GreenFunction.cpp
             src/RytovaKeldysh/RytovaKeldysh.cpp
    )
add_library(OBJECTS OBJECT ${_SOURCES} )
	


#adding all sources to target
link_libraries(OBJECTS)

#adding OpenMP
find_package(OpenMP)
link_libraries(OpenMP::OpenMP_CXX)

#adding MPI
if(EDUS_MPI)
  find_package(MPI REQUIRED)
  message("MPI include dir: ${MPI_CXX_INCLUDE_DIRS}")
  include_directories(${MPI_CXX_INCLUDE_DIRS})
  link_libraries(MPI::MPI_CXX)
  message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
  if(MPI_COMPILE_FLAGS)
    set(COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    set(LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif()

if(EDUS_HDF5)
    find_package(HDF5 REQUIRED)
    message(STATUS "HDF5 LIBRARY: ${HDF5_LIBRARIES}")
    message(STATUS "HDF5 parallel: ${HDF5_IS_PARALLEL}")
    message(STATUS "HDF5 include directories: ${HDF5_INCLUDE_DIRS}")
    include_directories(${HDF5_INCLUDE_DIRS})
    link_libraries(${HDF5_LIBRARIES})
endif()


#adding mkl (can it be replaced with lapack?)
set(MKL_ARCH intel64)

if(EDUS_MKL_THREAD)
  set(MKL_THREADING intel_thread)
else()
  set(MKL_THREADING sequential)
endif()

set(MKL_INTERFACE lp64)
find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
add_compile_options($<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
include_directories($<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
link_libraries($<LINK_ONLY:MKL::MKL>)
message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")


#fftw3

#find_package(fftw3 REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
#message(STATUS "PkgConfigFFTW -> ${FFTW_IMPORTED_TARGETS}")
include_directories(PkgConfig::FFTW)
link_libraries(PkgConfig::FFTW)
link_libraries(fftw3)
if(EDUS_FFTWTHREADS)
  link_libraries(fftw3_omp)
endif()
if (EDUS_MPI)
  link_libraries(fftw3_mpi)
endif()


# add the executable
add_executable(EDUS src/EDUS.cpp)

#speed executables
add_executable(Matrix_mult_speed speed-apps/Matrix_multiplication.cpp)
add_executable(BlockMatrix_mult_speed speed-apps/BlockMatrix_multiplication.cpp)
add_executable(FourierTransform_speed speed-apps/FourierTransform.cpp)

#Executables for test various classes
add_executable(Containers_test ci-test/test_mdarray.cpp) 
add_executable(Matrix_test ci-test/test_matrix.cpp)
add_executable(Eigen_test ci-test/test_eigen.cpp) 
add_executable(DE_test ci-test/test_DEsolver.cpp)
add_executable(DE2_test ci-test/test_DEsolver2.cpp)
#add_executable(fft-many_test ci-test/test_fft_many.cpp)
add_executable(fft_test ci-test/test_fft_mpi.cpp)
add_executable(Gradient_R_test ci-test/test_Gradient_R.cpp)
add_executable(NonEQ ci-test/test_NEQ_H0.cpp)
add_executable(NonEQ3 ci-test/test_NEQ_Grad.cpp)
add_executable(Bandstructure_test ci-test/test_bandstructure.cpp)
add_executable(Pseudoinv_test ci-test/test_pseudoinverse.cpp)

#build for ctest
enable_testing()
add_test(NAME Containers_test COMMAND Containers_test)
add_test(NAME Matrix_test COMMAND Matrix_test)
add_test(NAME Eigen_test COMMAND Eigen_test)
add_test(NAME DE_test COMMAND DE_test)
add_test(NAME DE2_test COMMAND DE2_test)
add_test(NAME Gradient_R_test COMMAND Gradient_R_test)
add_test(NAME Bandstructure_test COMMAND 
					bash -c "${CMAKE_CURRENT_BINARY_DIR}/Bandstructure_test ; \
					python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_CURRENT_BINARY_DIR}/BANDSTRUCTURE.txt ${CMAKE_SOURCE_DIR}/ci-test/Reference/BANDSTRUCTURE.txt")
add_test(NAME fft_test COMMAND fft_test)
add_test(NAME NonEQ_H0 COMMAND NonEQ)
add_test(NAME NonEQ_Grad COMMAND NonEQ3)

add_test(NAME hBN COMMAND
             bash -c "mkdir ${CMAKE_SOURCE_DIR}/CTEST; 
             rm -rf ${CMAKE_SOURCE_DIR}/CTEST/hBN/; \ 
             mkdir ${CMAKE_SOURCE_DIR}/CTEST/hBN/; \ 
             cd ${CMAKE_SOURCE_DIR}/CTEST/hBN/; \             
             ${CMAKE_CURRENT_BINARY_DIR}/EDUS  ${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN.json | tee output; \
             python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/hBN/Output/Population.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN/Population.txt
             python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/hBN/Output/Velocity.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN/Velocity.txt; \
             python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/hBN/Output/DM0.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN/DM0.txt")

add_test(NAME hBN_Coulomb COMMAND 
          bash -c "mkdir ${CMAKE_SOURCE_DIR}/CTEST; 
          rm -rf ${CMAKE_SOURCE_DIR}/CTEST/hBN_coulomb/; \ 
          mkdir ${CMAKE_SOURCE_DIR}/CTEST/hBN_coulomb/; \ 
          cd ${CMAKE_SOURCE_DIR}/CTEST/hBN_coulomb/; \             
          python3 ${PROJECT_SOURCE_DIR}/Postproces/RytovaKeldysh.py 40 40 1 ${PROJECT_SOURCE_DIR}/tb_models/hBN_gap7.25eV_a2.5A_tb.dat; \
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/hBN_coulomb/RytovaKeldysh.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN_coulomb/RytovaKeldysh.txt; \
          ${CMAKE_CURRENT_BINARY_DIR}/EDUS  ${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN_coulomb.json | tee output; \
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/hBN_coulomb/HF.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN_coulomb/HF.txt; \
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/hBN_coulomb/Output/Population.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN_coulomb/Population.txt; \
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/hBN_coulomb/Output/Velocity.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN_coulomb/Velocity.txt")

add_test(NAME MPI_hBN COMMAND
          bash -c "mkdir ${CMAKE_SOURCE_DIR}/CTEST; 
          rm -rf ${CMAKE_SOURCE_DIR}/CTEST/MPI_hBN/; \ 
          mkdir ${CMAKE_SOURCE_DIR}/CTEST/MPI_hBN/; \ 
          cd ${CMAKE_SOURCE_DIR}/CTEST/MPI_hBN/; \             
          ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS
          ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 10 ${CMAKE_CURRENT_BINARY_DIR}/EDUS  ${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN.json | tee output; \
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/MPI_hBN/Output/Population.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN/Population.txt
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/MPI_hBN/Output/Velocity.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN/Velocity.txt; \
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/MPI_hBN/Output/DM0.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN/DM0.txt")

add_test(NAME MPI_hBN_Coulomb COMMAND 
          bash -c "mkdir ${CMAKE_SOURCE_DIR}/CTEST; 
          rm -rf ${CMAKE_SOURCE_DIR}/CTEST/MPI_hBN_coulomb/; \ 
          mkdir ${CMAKE_SOURCE_DIR}/CTEST/MPI_hBN_coulomb/; \ 
          cd ${CMAKE_SOURCE_DIR}/CTEST/MPI_hBN_coulomb/; \             
          python3 ${PROJECT_SOURCE_DIR}/Postproces/RytovaKeldysh.py 40 40 1 ${PROJECT_SOURCE_DIR}/tb_models/hBN_gap7.25eV_a2.5A_tb.dat; \
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/MPI_hBN_coulomb/RytovaKeldysh.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN_coulomb/RytovaKeldysh.txt; \
          ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 10 ${CMAKE_CURRENT_BINARY_DIR}/EDUS ${PROJECT_SOURCE_DIR}/ci-test/inputs/hBN_coulomb.json | tee output; \
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/MPI_hBN_coulomb/HF.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN_coulomb/HF.txt; \
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/MPI_hBN_coulomb/Output/Population.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN_coulomb/Population.txt; \
          python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py ${CMAKE_SOURCE_DIR}/CTEST/MPI_hBN_coulomb/Output/Velocity.txt ${PROJECT_SOURCE_DIR}/ci-test/outputs/hBN_coulomb/Velocity.txt")
