#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <complex>

constexpr std::complex<double> im(0.,1.);
constexpr double pi = 3.14159265358979323846;


//CODATA values
constexpr double Bohr_value            = 5.2917721090380*1.e-11;//m
constexpr double ElectronVolt_value    = 1.602176634*1.e-19; //J
constexpr double Rydberg_value         = 2.1798723611035*1.e-18; //J
constexpr double AtomicUnitOfTime      = 2.4188843265857*1.e-17; //s
constexpr double Epsilon0              = 8.8541878128*1.e-12; //F m^-1
constexpr double E0                    = 5.14220674763*1.e+11;//V m^-1 (Atomic units of electric field) 
constexpr double FineStructure         = 7.2973525693*1.e-03;//pure units

constexpr double SpeedOfLight          = 1./FineStructure; //a.u.
constexpr double AtomicUnitOfIntensity = .5*SpeedOfLight*Epsilon0*E0*E0; //W/cm2

double c1,c2,c3;

double threshold = 1.e-08;

enum Space{k,R};
enum BandGauge{bloch, wannier};

#endif