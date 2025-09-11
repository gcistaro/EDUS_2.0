#include <map>
#include <string>
#include <cassert>
#include "ConvertUnits.hpp"


Unit::Unit(const Type& type_, const double& value_) : type(type_), value(value_){};

Unit Angstrom(LENGTH,1.e-10);
Unit NanoMeters(LENGTH, 1.e-09);
Unit ElectronVolt(ENERGY,ElectronVolt_value);
Unit Joule(ENERGY,1.);
Unit Wcm2(INTENSITY, 1.e+04);
Unit FemtoSeconds(TIME, 1.e-15);
Unit NanoSeconds(TIME, 1.e-9);
Unit MilliElectronVolt(ENERGY, ElectronVolt_value*1.e-3);
Unit Degrees(NullType, pi/180);

Unit AuIntensity(INTENSITY, AtomicUnitOfIntensity);
Unit AuTime(TIME, AtomicUnitOfTime);
Unit AuLength(LENGTH,Bohr_value);
Unit AuEnergy(ENERGY,2.*Rydberg_value);
Unit Rads(NullType, 1);
Unit NullUnit(NullType, 0.);

Unit unit( const std::string& to_unit)
{
    if( to_unit == "angstrom"    )    {return Angstrom;           }
    if( to_unit == "nanometers"  )    {return NanoMeters;         }
    if( to_unit == "nm"          )    {return NanoMeters;         }
    if( to_unit == "electronvolt")    {return ElectronVolt;       }
    if( to_unit == "joule"       )    {return Joule;              }
    if( to_unit == "wcm2"        )    {return Wcm2;               }
    if( to_unit == "femtoseconds")    {return FemtoSeconds;       }
    if( to_unit == "fs"          )    {return FemtoSeconds;       }
    if( to_unit == "auintensity" )    {return AuIntensity;        }
    if( to_unit == "autime"      )    {return AuTime;             }
    if( to_unit == "aulength"    )    {return AuLength;           }
    if( to_unit == "auenergy"    )    {return AuEnergy;           }
    if( to_unit == "ns"          )    {return NanoSeconds;        }
    if( to_unit == "meV"         )    {return MilliElectronVolt;  }
    if( to_unit == "degrees"     )    {return Degrees;            }
    if( to_unit == "rads"        )    {return Rads;               }
    return NullUnit;
}
