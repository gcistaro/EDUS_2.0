#include "mdarray.h"

int main()
{
    //void constructor
    mdarray<double,2> VoidConstructed;

    //move assignment
    auto MoveAssigned = mdarray<double,2>({5,5});

    //copy assignemnt
    auto CopyAssigned = MoveAssigned;

    //copy constructor
    auto CopyConstructed = mdarray<double,2>(MoveAssigned);

    //access const element
    auto i = MoveAssigned[0];
    i = MoveAssigned(0,1);

    //access nonconst element
    MoveAssigned[0] = i;
    MoveAssigned(0,2) = i;
}
