#include "GlobalFunctions.hpp"

std::string LatticeVectors(const Space& space) 
{ 
    return ( space == k ? "k_LatticeVectors" : "R_LatticeVectors" ); 
};


bool
file_exists(std::string file_name)
{
    std::ifstream ifs(file_name.c_str());
    return ifs.is_open();
}


void dump_json_in_h5( const nlohmann::json& data__, const std::string& name__ )
{
#ifdef EDUS_HDF5
    HDF5_tree fout(name__, hdf5_access_t::read_write);
    fout.create_node("json");

    fout["json"].write("tb_file", data__["tb_file"].template get<std::string>());
    //fout["json"].write("fermienergy", 
    //        {data__["fermienergy"][0].template get<double>(), data__["fermienergy"][1].template get<std::string>()});
    fout["json"].write("coulomb", data__["coulomb"].template get<int>());
    auto grid = std::vector<int>(3);
    grid = {data__["grid"][0].template get<int>(), 
            data__["grid"][1].template get<int>(), 
            data__["grid"][2].template get<int>()}; 
    fout["json"].write("grid", grid);
    fout["json"].write("filledbands", data__["filledbands"].template get<int>());
    fout["json"].write("dt", data__["dt"][0].template get<double>());
    fout["json"].write("printresolution", data__["printresolution"].template get<int>());
    fout["json"].write("initialtime", data__["initialtime"][0].template get<double>());
    fout["json"].write("finaltime", data__["finaltime"][0].template get<double>());
    fout["json"].write("solver", data__["solver"].template get<std::string>());
    for ( int ilaser = 0; ilaser < int( data__["lasers"].size() ); ++ilaser ) {
        auto lasername = "laser"+std::to_string(ilaser);
        fout["json"].create_node(lasername);
        fout["json"][lasername].write("intensity", data__["lasers"][ilaser]["intensity"][0].template get<double>());
        fout["json"][lasername].write("frequency", data__["lasers"][ilaser]["frequency"][0].template get<double>());
        fout["json"][lasername].write("wavelength", data__["lasers"][ilaser]["wavelength"][0].template get<double>());
        fout["json"][lasername].write("cycles", data__["lasers"][ilaser]["cycles"].template get<int>());
        auto grid = std::vector<int>(3);
        grid = { data__["lasers"][ilaser]["polarization"][0].template get<int>(), 
                 data__["lasers"][ilaser]["polarization"][1].template get<int>(), 
                 data__["lasers"][ilaser]["polarization"][2].template get<int>()}; 

        fout["json"][lasername].write("polarization", data__["lasers"][ilaser]["cycles"].template get<int>());
    }
#endif
}