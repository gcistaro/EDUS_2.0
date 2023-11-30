#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <complex>

constexpr std::complex<double> im(0.,1.);
constexpr double pi = 3.14159265358979323846;


//CODATA values
constexpr double Bohr_value         = 5.2917721090380*1.e-11;//m
constexpr double ElectronVolt_value = 1.602176634*1.e-19; //J
constexpr double Rydberg_value      = 2.1798723611035*1.e-18; //J
constexpr double AtomicUnitOfTime   = 2.4188843265857*1.e-17; //s

constexpr double SpeedOfLight       = 299792458/Bohr_value*AtomicUnitOfTime;

double c1,c2,c3;

double threshold = 1.e-08;

#endif