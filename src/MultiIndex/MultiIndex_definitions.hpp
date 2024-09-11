#include "MultiIndex.hpp"
#include <iostream>

template<size_t dim> 
void MultiIndex<dim>::TotalSizeAndOffset()
{
    TotalSize = 1;
    for(size_t idim=0; idim<dim; idim++){
        TotalSize *= Size[idim]; 
    }
    Offset[dim-1] = 1;
    for(int idim=int(dim)-2; idim>=0; idim--){
        Offset[idim] = Size[idim+1]*Offset[idim+1];
    }
}

template<size_t dim> 
MultiIndex<dim>::MultiIndex(const std::array<int,dim>& Size__)
{
    initialize(Size__);
}

template<size_t dim> 
void MultiIndex<dim>::initialize(const std::array<int,dim>& Size__) 
{
    Size = Size__;
    this->TotalSizeAndOffset();
}

template<size_t dim> 
template <typename... Args>
inline int MultiIndex<dim>::oneDindex(Args... args)
{
    return const_cast<const MultiIndex<dim>&>(*this).oneDindex(args...);
}

template<size_t dim> 
template <typename... Args>
inline int MultiIndex<dim>::oneDindex(Args... args) const
{
    assert( dim == sizeof...(args) );
    std::array<int,dim> i = {args...};

    int oneDindex = i[dim-1]*Offset[dim-1];

    for(int idim=int(dim)-2; idim>=0; idim--){
        oneDindex += Offset[idim]*i[idim];               
    }
    return oneDindex;
}


template<size_t dim> 
inline std::array<int,dim> MultiIndex<dim>::nDindex(const int& index)
{
    return const_cast<const MultiIndex<dim>&>(*this).nDindex(index);
}

template<size_t dim> 
inline std::array<int,dim> MultiIndex<dim>::nDindex(const int& index) const
{
    assert(index >= 0 && index < TotalSize);
    std::array<int,dim> nDindex;
    int remainder = index;
    for(int idim=0; idim<dim; idim++){
        nDindex[idim] = remainder/Offset[idim];
        remainder = remainder%Offset[idim];
    }
    return nDindex;
}

