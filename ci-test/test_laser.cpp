#include "Laser/Laser.hpp"

int main()
{
    std::cout << "Testing laser...\n";
    Laser laser;
    laser.set_InitialTime(0., FemtoSeconds);
    laser.set_Intensity(1.e+05, Wcm2);
    laser.set_Polarization(Coordinate({1,0,0}));
    laser.set_Lambda(800, NanoMeters);
    laser.set_NumberOfCycles(10);

    std::ofstream os_Laser("Laser_test.txt");
    //os_Laser << laser;
    os_Laser.close();
}
