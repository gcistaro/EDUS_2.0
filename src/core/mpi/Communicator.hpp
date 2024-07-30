//this file is adapted from sirius code
#ifndef COMMUNICATORS_HPP
#define COMMUNICATORS_HPP


#ifdef NEGF_MPI

#include "mpi.h"

namespace mpi
{

class Communicator {
    private:
        MPI_Comm communicator_{MPI_COMM_NULL};
        int rank_{-1};
        int size_{-1};

        void init()
        {
            assert(communicator_ != MPI_COMM_NULL);
            MPI_Comm_rank(communicator_, &rank_);
            MPI_Comm_size(communicator_, &size_);
        };
    
    public:
        //Default
        Communicator(){};

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
    
        /// MPI global initialization.
        static void initialize(int required__)
        {
            int provided;

            MPI_Init_thread(NULL, NULL, required__, &provided);
            MPI_Query_thread(&provided);
            if ((provided < required__) && (Communicator::world().rank() == 0)) {
                std::printf("Warning! Required level of thread support is not provided.\n");
                std::printf("provided: %d \nrequired: %d\n", provided, required__);
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