//this function must be called at the beginning of any miniapp to allow for 
//initialization of NEGF important quantities, such as:
//- header printing with git infos
//- MPI_COMM_WORLD
//- fftw initialization
#ifndef INITIALIZE_HPP
#define INITIALIZE_HPP

#include <sstream>
#include <fstream>

#include "core/print_header.hpp"
#include "core/mpi/Communicator.hpp"
#include "fftw3-mpi.h"


#ifdef NEGF_MPI
extern mpi::Communicator kpool_comm;
extern mpi::Communicator band_comm;
extern int NumberKpools;
#endif

#define variable(x)  (#x)

void initialize();

#endif