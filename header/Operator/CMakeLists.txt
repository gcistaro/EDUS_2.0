cmake_minimum_required(VERSION 3.16) #last feature: target_precompile_header

# set the project name and version
project(NEGF VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


#add flags
SET(CXX_COVERAGE_COMPILE_FLAGS "-g -fpermissive")
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS})


# add the executable
add_executable(NEGF main.cpp)

#mkl
find_package(MKL CONFIG REQUIRED)
target_compile_options(NEGF PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(NEGF PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(NEGF PUBLIC $<LINK_ONLY:MKL::MKL>)

#headers
###target_include_directories(NEGF PUBLIC header/)
