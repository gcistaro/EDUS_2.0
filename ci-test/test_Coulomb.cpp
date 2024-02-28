#include "Model/Model.hpp"
#include "Coulomb/Coulomb.hpp"

int main()
{
    Material material("/home/gcistaro/NEGF/tb_models/non_dummy_gap");
    Coulomb coulomb(material.r, 2);

    for(int i=0; i<1000; i++){
        auto r = Coordinate<R>(double(i)/1000., 0., 0.);
        std::cout << double(i)/100.*10 << " " << (coulomb.W(r)).real() << std::endl;
    }

}