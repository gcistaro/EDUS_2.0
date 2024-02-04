#include "Operator.hpp"

int main()
{
    size_t nb = 7;
    BandIndex bi(nb);
    for(int i=0; i<nb; ++i){
        std::cout << i << " " << bi.StartingIndex(i)<<std::endl;
    }

    std::cout << "bi.get_oneDNumberOfBands(): " << bi.get_oneDNumberOfBands() << std::endl;
    for(int i=0; i<bi.get_oneDNumberOfBands(); ++i){
        auto index = bi.twoDband(i);
        std::cout << i << " " << index.first << " " << index.second << std::endl;
    }

}