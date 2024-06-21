cmake_minimum_required(VERSION 3.16) #last feature: target_precompile_header

# set the project name and version
project(NEGF VERSION 1.0)

option(NEGF_PROFILE "Profiling function timings" ON)

set(CMAKE_VERBOSE_MAKEFILE ON)
# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /lib/x86_64-linux-gnu/)
#set(CMAKE_REQUIRED_FLAGS -fconcepts)
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
set(GCC_FLAG "-fconcepts -g")
add_definitions(${GCC_FLAG} -DNEGF_PROFILE=ON)


#git version
# check if git command exists
find_program(GIT_EXE NAMES git)

# generate version header
string(TIMESTAMP NEGF_TIMESTAMP "%Y-%m-%d %H:%M:%S")
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

#compile sources as library INTERFACE (not a real library, just objects to be reused)
include_directories("${PROJECT_SOURCE_DIR}/src/")

set(_SOURCES src/GlobalFunctions.cpp
             src/ConvertUnits.cpp
	     src/ostream.cpp
             src/StreamFile.cpp
             src/Geometry/Matrix.cpp 
             src/Geometry/Coordinate.cpp
	     src/MeshGrid/MeshGrid.cpp
             src/Wannier/Wannier.cpp
             src/Laser/Laser.cpp
             src/fftPair/fftPair.cpp
	     src/kGradient/kGradient.cpp
             src/Simulation/Simulation.cpp
             src/core/profiler.cpp 
             src/core/print_timing.cpp 
             src/core/rt_graph.cpp
			 )

find_package(OpenMP)
link_libraries(OpenMP::OpenMP_CXX)
                      
add_library(OBJECTS OBJECT ${_SOURCES})
#target_sources(OBJECTS PUBLIC ${_SOURCES})
link_libraries(OBJECTS)

#mkl
set(MKL_ARCH intel64)
set(MKL_THREADING sequential)
set(MKL_INTERFACE lp64)
find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
add_compile_options($<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
include_directories($<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
link_libraries($<LINK_ONLY:MKL::MKL>)

message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")

# add the executable
add_executable(NEGF main.cpp)


#headers
target_precompile_headers(NEGF PUBLIC 
	                 src/Geometry/Matrix.hpp
	                 src/Geometry/Vector.hpp
 			   src/mdContainers/mdContainers.hpp
 			 )

#testing...
#mdcontainers
add_executable(Containers_test ci-test/test_mdarray.cpp) 

#Matrix
add_executable(Matrix_test ci-test/test_matrix.cpp)

#Laser
add_executable(Laser_test ci-test/test_laser.cpp)

#Eigenvalues
add_executable(Eigen_test ci-test/test_eigen.cpp) 

#RK
add_executable(RK_test ci-test/test_RK.cpp)

#RK2
add_executable(RK2_test ci-test/test_RK2.cpp)

#fft
add_executable(fft_test ci-test/test_fft.cpp)
#find_package(fftw3 REQUIRED)
target_link_libraries(fft_test PUBLIC fftw3)

#fft-many
add_executable(fft-many_test ci-test/test_fft_many.cpp)
#find_package(fftw3 REQUIRED)
target_link_libraries(fft-many_test PUBLIC fftw3)

#Gradient
add_executable(Gradient_test ci-test/test_Gradient.cpp)
target_link_libraries(Gradient_test PUBLIC fftw3)

#fft+gradient
add_executable(fft_Gradient_test ci-test/test_fft+gradient.cpp)
target_link_libraries(fft_Gradient_test PUBLIC fftw3)

#TrivialBands 
#add_executable(TrivialH_test ci-test/test_TrivialHamiltonian.cpp)
#find_package(fftw3 REQUIRED)
#target_link_libraries(TrivialH_test PUBLIC fftw3)

#Noneq propagation 
add_executable(NonEQ ci-test/test_NEQpop2.cpp)
#find_package(fftw3 REQUIRED)
target_link_libraries(NonEQ PUBLIC fftw3)

#simulation 
add_executable(simulation_test ci-test/test_simulation.cpp)
#find_package(fftw3 REQUIRED)
target_link_libraries(simulation_test PUBLIC fftw3)

#Coulomb 
#add_executable(Coulomb_test ci-test/test_Coulomb.cpp)
#find_package(fftw3 REQUIRED)
#target_link_libraries(Coulomb_test PUBLIC fftw3)

#bandstructure
add_executable(Bandstructure_test ci-test/test_bandstructure.cpp)
#find_package(fftw3 REQUIRED)
target_link_libraries(Bandstructure_test PUBLIC fftw3)

#pseudoinverse
add_executable(Pseudoinv_test ci-test/test_pseudoinverse.cpp)
#find_package(fftw3 REQUIRED)
target_link_libraries(Pseudoinv_test PUBLIC fftw3)

#build for ctest
enable_testing()
add_test(NAME Containers_test COMMAND Containers_test)
add_test(NAME Matrix_test COMMAND Matrix_test)
add_test(NAME Laser_test COMMAND Laser_test)
add_test(NAME Eigen_test COMMAND Eigen_test)
add_test(NAME RK_test COMMAND RK_test)
add_test(NAME RK2_test COMMAND RK2_test)
add_test(NAME fft_test COMMAND fft_test)
add_test(NAME fft-many_test COMMAND fft-many_test)
add_test(NAME fft_Gradient_test COMMAND fft_Gradient_test)
#add_test(NAME TrivialH_test COMMAND TrivialH_test)
add_test(NAME Bandstructure_test COMMAND 
					bash -c "${CMAKE_CURRENT_BINARY_DIR}/Bandstructure_test ; \
					python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py BANDSTRUCTURE.txt ${CMAKE_SOURCE_DIR}/ci-test/Reference/BANDSTRUCTURE.txt")
add_test(NAME simulation_test COMMAND simulation_test)
add_test(NAME NonEQ_test COMMAND NonEQ_test)
add_test(NAME Coulomb_test COMMAND Coulomb_test)

