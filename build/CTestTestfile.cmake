# CMake generated Testfile for 
# Source directory: /home/miguelsa/Desktop/madrid_phd/NEGF
# Build directory: /home/miguelsa/Desktop/madrid_phd/NEGF/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(Containers_test "/home/miguelsa/Desktop/madrid_phd/NEGF/build/Containers_test")
set_tests_properties(Containers_test PROPERTIES  _BACKTRACE_TRIPLES "/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;233;add_test;/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;0;")
add_test(MultiIndex_test "/home/miguelsa/Desktop/madrid_phd/NEGF/build/MultiIndex_test")
set_tests_properties(MultiIndex_test PROPERTIES  _BACKTRACE_TRIPLES "/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;234;add_test;/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;0;")
add_test(Matrix_test "/home/miguelsa/Desktop/madrid_phd/NEGF/build/Matrix_test")
set_tests_properties(Matrix_test PROPERTIES  _BACKTRACE_TRIPLES "/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;235;add_test;/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;0;")
add_test(Eigen_test "/home/miguelsa/Desktop/madrid_phd/NEGF/build/Eigen_test")
set_tests_properties(Eigen_test PROPERTIES  _BACKTRACE_TRIPLES "/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;237;add_test;/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;0;")
add_test(RK_test "/home/miguelsa/Desktop/madrid_phd/NEGF/build/RK_test")
set_tests_properties(RK_test PROPERTIES  _BACKTRACE_TRIPLES "/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;238;add_test;/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;0;")
add_test(RK2_test "/home/miguelsa/Desktop/madrid_phd/NEGF/build/RK2_test")
set_tests_properties(RK2_test PROPERTIES  _BACKTRACE_TRIPLES "/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;239;add_test;/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;0;")
add_test(Gradient_R_test "/home/miguelsa/Desktop/madrid_phd/NEGF/build/Gradient_R_test")
set_tests_properties(Gradient_R_test PROPERTIES  _BACKTRACE_TRIPLES "/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;242;add_test;/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;0;")
add_test(Bandstructure_test "bash" "-c" "/home/miguelsa/Desktop/madrid_phd/NEGF/build/Bandstructure_test ; 					python3 /home/miguelsa/Desktop/madrid_phd/NEGF/ci-test/compare.py BANDSTRUCTURE.txt /home/miguelsa/Desktop/madrid_phd/NEGF/ci-test/Reference/BANDSTRUCTURE.txt")
set_tests_properties(Bandstructure_test PROPERTIES  _BACKTRACE_TRIPLES "/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;244;add_test;/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;0;")
add_test(fft_mpi_test "/home/miguelsa/Desktop/madrid_phd/NEGF/build/fft_mpi_test")
set_tests_properties(fft_mpi_test PROPERTIES  _BACKTRACE_TRIPLES "/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;247;add_test;/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;0;")
add_test(NonEQ "/home/miguelsa/Desktop/madrid_phd/NEGF/build/NonEQ")
set_tests_properties(NonEQ PROPERTIES  _BACKTRACE_TRIPLES "/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;249;add_test;/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;0;")
add_test(NonEQ3 "/home/miguelsa/Desktop/madrid_phd/NEGF/build/NonEQ3")
set_tests_properties(NonEQ3 PROPERTIES  _BACKTRACE_TRIPLES "/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;250;add_test;/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;0;")
add_test(RytKel "/home/miguelsa/Desktop/madrid_phd/NEGF/build/RytKel")
set_tests_properties(RytKel PROPERTIES  _BACKTRACE_TRIPLES "/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;251;add_test;/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;0;")
add_test(Overhead "/home/miguelsa/Desktop/madrid_phd/NEGF/build/Overhead")
set_tests_properties(Overhead PROPERTIES  _BACKTRACE_TRIPLES "/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;252;add_test;/home/miguelsa/Desktop/madrid_phd/NEGF/CMakeLists.txt;0;")
