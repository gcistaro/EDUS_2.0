cmake_minimum_required(VERSION 3.16) #last feature: target_precompile_header

# set the project name and version
project(NEGF VERSION 1.0)


option(NEGF_PROFILE "Profiling function timings" ON)

set(CMAKE_VERBOSE_MAKEFILE ON)
# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(CMAKE_REQUIRED_FLAGS -fconcepts)

set(GCC_FLAG "-fconcepts -g")
add_definitions(${GCC_FLAG} -DNEGF_PROFILE=ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

#compile sources as library INTERFACE (not a real library, just objects to be reused)
set(_SOURCES header/core/profiler.cpp 
             header/Geometry/Matrix.cpp 
             header/core/print_timing.cpp 
             header/core/rt_graph.cpp)
add_library(OBJECTS INTERFACE)
target_sources(OBJECTS PUBLIC ${_SOURCES})


set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /lib/x86_64-linux-gnu/)
# add the executable
add_executable(NEGF main.cpp)
target_link_libraries(NEGF PUBLIC OBJECTS)

#mkl
set(MKL_ARCH intel64)
set(MKL_THREADING sequential)
set(MKL_INTERFACE lp64)
find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
target_compile_options(NEGF PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(NEGF PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(NEGF PUBLIC $<LINK_ONLY:MKL::MKL>)
target_include_directories(NEGF PUBLIC header/)

message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")

#headers
target_precompile_headers(NEGF PUBLIC 
	                 header/Geometry/Matrix.hpp
	                 header/Geometry/Vector.hpp
 			 header/mdContainers/mdContainers.hpp
 			 )

#testing...
#mdcontainers
add_executable(Containers_test ci-test/test_mdarray.cpp) #${_SOURCES})
target_include_directories(Containers_test PUBLIC header/)

#Matrix
add_executable(Matrix_test ci-test/test_matrix.cpp)
target_compile_options(Matrix_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(Matrix_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(Matrix_test PUBLIC $<LINK_ONLY:MKL::MKL>)
target_include_directories(Matrix_test PUBLIC header/)
target_link_libraries(Matrix_test PUBLIC OBJECTS)

#Eigenvalues
add_executable(Eigen_test ci-test/test_eigen.cpp) 
target_compile_options(Eigen_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(Eigen_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(Eigen_test PUBLIC $<LINK_ONLY:MKL::MKL>)
target_include_directories(Eigen_test PUBLIC header/)
target_link_libraries(Eigen_test PUBLIC OBJECTS)


#RK
add_executable(RK_test ci-test/test_RK.cpp)
target_include_directories(RK_test PUBLIC header/)

#fft
add_executable(fft_test ci-test/test_fft.cpp)
#find_package(fftw3 REQUIRED)
target_link_libraries(fft_test fftw3)
target_include_directories(fft_test PUBLIC header/)

#simulation 
add_executable(simulation_test ci-test/test_simulation.cpp)
#find_package(fftw3 REQUIRED)
target_link_libraries(simulation_test PUBLIC fftw3)
target_compile_options(simulation_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(simulation_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(simulation_test PUBLIC $<LINK_ONLY:MKL::MKL>)
target_include_directories(simulation_test PUBLIC header/)
target_link_libraries(simulation_test PUBLIC OBJECTS)

#simulation 
add_executable(Coulomb_test ci-test/test_Coulomb.cpp)
#find_package(fftw3 REQUIRED)
target_link_libraries(Coulomb_test PUBLIC fftw3)
target_compile_options(Coulomb_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(Coulomb_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(Coulomb_test PUBLIC $<LINK_ONLY:MKL::MKL>)
target_include_directories(Coulomb_test PUBLIC header/)
target_link_libraries(Coulomb_test PUBLIC OBJECTS)

#bandstructure
add_executable(Bandstructure_test ci-test/test_bandstructure.cpp)
#find_package(fftw3 REQUIRED)
target_link_libraries(Bandstructure_test PUBLIC fftw3)
target_compile_options(Bandstructure_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(Bandstructure_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(Bandstructure_test PUBLIC $<LINK_ONLY:MKL::MKL>)
target_include_directories(Bandstructure_test PUBLIC header/)
target_link_libraries(Bandstructure_test PUBLIC OBJECTS)

#build for ctest
enable_testing()
add_test(NAME Containers_test COMMAND Containers_test)
add_test(NAME Matrix_test COMMAND Matrix_test)
add_test(NAME Eigen_test COMMAND Eigen_test)
add_test(NAME RK_test COMMAND RK_test)
add_test(NAME fft_test COMMAND fft_test)
add_test(NAME simulation_test COMMAND simulation_test)
add_test(NAME Coulomb_test COMMAND Coulomb_test)
add_test(NAME Bandstructure_test COMMAND 
					bash -c "${CMAKE_CURRENT_BINARY_DIR}/Bandstructure_test ; \
					python3 ${CMAKE_SOURCE_DIR}/ci-test/compare.py BANDSTRUCTURE.txt ${CMAKE_SOURCE_DIR}/ci-test/Reference/BANDSTRUCTURE.txt")
