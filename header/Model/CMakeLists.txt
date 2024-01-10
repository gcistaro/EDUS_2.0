cmake_minimum_required(VERSION 3.16) #last feature: target_precompile_header

# set the project name and version
project(NEGF VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# add the executable
add_executable(NEGF bandstructure.cpp)

#mkl
find_package(MKL CONFIG REQUIRED)
target_include_directories(NEGF PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(NEGF PUBLIC $<LINK_ONLY:MKL::MKL>)

#headers
target_include_directories(NEGF PUBLIC ../header/Model
/home/gcistaro/NEGF/header/fftPair
/home/gcistaro/NEGF/header/Geometry
/home/gcistaro/NEGF/header/Laser
/home/gcistaro/NEGF/header/LinearAlgebra
/home/gcistaro/NEGF/header/mdContainers
/home/gcistaro/NEGF/header/MeshGrid
/home/gcistaro/NEGF/header/Model
/home/gcistaro/NEGF/header/Operator
/home/gcistaro/NEGF/header/RungeKutta
/home/gcistaro/NEGF/header/Simulation
/home/gcistaro/NEGF/header/TimeGrid
/home/gcistaro/NEGF/header/Wannier
/home/gcistaro/NEGF/header/)


#target_precompile_headers(NEGF PUBLIC 
#	                 header/Matrix.hpp
#	                 header/Geometry.hpp
#			 header/MeshGrid.hpp
#			 header/Tensor.hpp)
