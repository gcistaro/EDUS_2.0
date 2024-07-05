
#ifndef MULTIINDEX_HPP
#define MULTIINDEX_HPP
#include <array>
#include <vector>
#include <cassert>
/*
   This class can be used to get oneD indices of multiarrays
   or to get nD indices
*/
template<size_t dim> 
class MultiIndex 
{
    private:
        std::array<size_t,dim> Size;
        std::array<size_t,dim> Offset;
        int TotalSize;
    public:
        MultiIndex(){};
        MultiIndex(const MultiIndex& ToBeCopied) = default;
        MultiIndex& operator=(const MultiIndex& ToBeCopied) = default;
        
        MultiIndex(MultiIndex&& ToBeMoved) = default;
        MultiIndex& operator=(MultiIndex&& ToBeMoved) = default;

        MultiIndex(const std::array<size_t,dim>& Dimension__);        
        void initialize(const std::array<size_t,dim>& Dimension__);   
        void TotalSizeAndOffset();
       
        template <typename... Args>
        inline int oneDindex(Args... args);

        inline std::array<int,dim> nDindex(const int& index);
};

#include "MultiIndex_definitions.hpp"
#endif

