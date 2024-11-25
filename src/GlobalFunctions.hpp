#include "Constants.hpp"
#include <fstream>

#ifdef EDUS_HDF5
#include "core/hdf5/hdf5_tree.hpp"
#endif
#include <string>
#include "Json/json.hpp"

std::string LatticeVectors(const Space& space);
bool file_exists(std::string file_name);
void dump_json_in_h5( const nlohmann::json& data__, const std::string& name__ );
