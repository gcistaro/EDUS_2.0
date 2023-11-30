#include "Laser.hpp"

int main()
{
    TimeGrid timegrid(0., 100., 0.1);
    //std::cout << timegrid;
    Laser laser;
    laser.set_Omega(1.);
    laser.set_InitialTime(0.);
    laser.set_NumberOfCycles(10);
    laser.set_Amplitude(110.);
    laser.set_Polarization(Vector<k>(1.,0.,0.,"Cartesian"));
    std::cout << laser(0.1);
}