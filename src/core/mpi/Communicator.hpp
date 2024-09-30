//this file is adapted from sirius code
#ifndef COMMUNICATORS_HPP
#define COMMUNICATORS_HPP


#ifdef NEGF_MPI

#include <complex>
#include "mpi.h"
#include <cassert>
#include <iostream>

namespace mpi
{

int tag ( const int& sender, const int& receiver );

class Communicator {
    private:
        MPI_Comm communicator_{MPI_COMM_NULL};
        int rank_{-1};
        int size_{-1};
    
    public:
        //Default
        Communicator(){};

        void init()
        {
            assert(communicator_ != MPI_COMM_NULL);
            MPI_Comm_rank(communicator_, &rank_);
            MPI_Comm_size(communicator_, &size_);
        };

        //From raw communicator
        Communicator(MPI_Comm mpi_comm__)
            : communicator_(mpi_comm__)
        {
            init();
        };

        /// Rank of actual process inside the communicator
        inline int rank() const
        {
            return rank_;
        };
    
        /// Size of actual communicator
        inline int size() const
        {
            return size_;
        };

        //communicator of the object
        inline MPI_Comm&
        communicator() 
        {
            return communicator_;
        }

        inline const MPI_Comm&
        communicator() const
        {
            return communicator_;
        }

        //split from a communicator with given color
        void generate( const Communicator& father_comm, int color )
        {
            MPI_Comm_split(father_comm.communicator(), color, father_comm.rank(), &communicator_);
            init();
        }

        //send message -- to be optimized with any type 
        void send( const int& message, const int& receiver ) const
        {
            MPI_Send( &message, 1, MPI_INT, receiver, mpi::tag(rank(), receiver), communicator_);
        }

        //receive message -- to be optimized with any type 
        void receive( int& message, const int& sender ) const
        {
            MPI_Recv( &message, 1, MPI_INT, sender, mpi::tag(sender, rank()), communicator_, MPI_STATUS_IGNORE);
        }

        //send message -- to be optimized with any type 
        void isend( const int& message, const int& receiver ) const
        {
            MPI_Request request;
            MPI_Isend( &message, 1, MPI_INT, receiver, mpi::tag(rank(), receiver), communicator_, &request);
            MPI_Wait(&request, MPI_STATUS_IGNORE);
        }

        //receive message -- to be optimized with any type 
        void ireceive( int& message, const int& sender ) const
        {
            MPI_Request request;
            MPI_Irecv( &message, 1, MPI_INT, sender, mpi::tag(sender, rank()), communicator_, &request);
        }

        //gather
        void reduce( std::complex<double>* to_send, std::complex<double>* to_recv, int count, const MPI_Op& mpi_op, int root)
        {
            MPI_Reduce( to_send, to_recv, count, MPI_C_DOUBLE_COMPLEX, mpi_op, root, communicator_);
        }

        /// MPI global initialization.
        static void initialize(int required__)
        {
            int provided;

            MPI_Init_thread(NULL, NULL, required__, &provided);
            MPI_Query_thread(&provided);
            if ((provided < required__) && (Communicator::world().rank() == 0)) {
                printf("Warning! Required level of thread support is not provided.\n");
                printf("provided: %d \nrequired: %d\n", provided, required__);
            }
        };

        // MPI global finalization
        static void finalize()
        {
            MPI_Finalize();
        };
        //retrieve MPI_COMM_WORLD communicator
        static const Communicator& world()
        {
            static Communicator comm(MPI_COMM_WORLD);
            return comm;
        };
};


}//end namespace mpi
#endif//NEGF_MPI
#endif//COMMUNICATORS_HPP