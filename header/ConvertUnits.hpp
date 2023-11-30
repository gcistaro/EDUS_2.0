#ifndef CONVERTUNITS_HPP
#define CONVERTUNITS_HPP

#include <map>
#include <string>
#include <cassert>

#include "Constants.hpp"


enum Type{LENGTH, TIME, ENERGY};

struct Unit{
    Type type;
    double value;
    
    Unit(const Type& type_, const double& value_) : type(type_), value(value_){};
};

Unit Angstrom(LENGTH,1.e-10);
Unit Bohr(LENGTH,Bohr_value);
Unit Rydberg(ENERGY,Rydberg_value);
Unit ElectronVolt(ENERGY,ElectronVolt_value);
Unit Joule(ENERGY,1.);


template<typename T>
T Convert(const T& ConvertableValue, const Unit& InputUnit, const Unit& OutputUnit)
{
    assert(InputUnit.type == OutputUnit.type);
    return ConvertableValue/InputUnit.value*OutputUnit.value; 
}

//we need a concept with iterators!!
template<typename T>
void Convert_iterable(T& ConvertableTensor, const Unit& InputUnit, const Unit& OutputUnit)
{
    for(auto& ToConvert : ConvertableTensor){
        ToConvert = Convert(ToConvert, InputUnit, OutputUnit);
    }
}



#endif