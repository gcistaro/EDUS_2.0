cmake_minimum_required(VERSION 3.16) #last feature: target_precompile_header

# set the project name and version
project(NEGF VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# add the executable
add_executable(NEGF main.cpp)

#mkl
set(MKL_ARCH intel64)
set(MKL_THREADING sequential)
set(MKL_INTERFACE lp64)
find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
target_compile_options(NEGF PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(NEGF PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(NEGF PUBLIC $<LINK_ONLY:MKL::MKL>)
message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")

#headers
target_include_directories(NEGF PUBLIC ./header/)
target_precompile_headers(NEGF PUBLIC 
	                 Matrix.hpp
	                 Vector.hpp
			 )
