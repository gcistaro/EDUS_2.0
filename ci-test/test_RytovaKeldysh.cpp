#include <cassert>
#include "Model/Model.hpp"
#include "initialize.hpp"
#include "RytovaKeldysh/RytovaKeldysh.hpp"




int main()
{   
    initialize();
    //------------------------------------------initialize obj RytovaKeldysh------------------------------------------------------
    auto material = Material("/home/gcistaro/NEGF/tb_models/hBN_gap7.25eV_a2.5A");
    auto MasterRgrid = std::make_shared<MeshGrid>(get_GammaCentered_grid(MeshGrid(R, {4,4,1})));
    RytovaKeldysh RytKeld(material.r, 2, MasterRgrid);
    
    std::ofstream os("RytovaKeldysh.txt");
    for(int i=0; i<100; i++) {
        os << i/10. << "  " << RytKeld.Potential(Coordinate(i/10.,0.,0.)).real()<< std::endl;
    }
    os.close();


    //------------------------------------------get wanted initial condition------------------------------------------------
}
