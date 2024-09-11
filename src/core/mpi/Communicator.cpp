#include "Communicator.hpp"

namespace mpi {
    int tag ( const int& sender, const int& receiver )
    {
        //static int counter = 0;
        //counter++;
        return sender + 2000*receiver; //+ counter*10000;
    }

}