#include <map>
#include <string>
#include <cassert>

#include "ConvertUnits.hpp"


Unit::Unit(const Type& type_, const double& value_) : type(type_), value(value_){};

Unit Angstrom(LENGTH,1.e-10);
Unit NanoMeters(LENGTH, 1.e-09);
Unit ElectronVolt(ENERGY,ElectronVolt_value);
Unit Joule(ENERGY,1.);
Unit Wcm2(INTENSITY, 1.);
Unit FemtoSeconds(TIME, 1.e-15);

Unit AuIntensity(INTENSITY, AtomicUnitOfIntensity);
Unit AuTime(TIME, AtomicUnitOfTime);
Unit AuLength(LENGTH,Bohr_value);
Unit AuEnergy(ENERGY,2.*Rydberg_value);
