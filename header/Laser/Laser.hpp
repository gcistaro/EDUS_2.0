#ifndef LASER_HPP
#define LASER_HPP

//#include <math.h>

#include "Constants.hpp"
#include "ConvertUnits.hpp"
#include "Geometry/Coordinate.hpp"
#include "TimeGrid/TimeGrid.hpp"

struct period
{   
    double value;
    period(const double& value_) : value(value_){};
};

struct frequency
{
    double value;
    frequency(const double& value_) : value(value_){};
};

struct wavelength
{
    double value;
    wavelength(const double& value_) : value(value_){};
};
//next three are the equations we will use for the others:
double get_period(const frequency& omega_) 
{
    return 2.*pi/omega_.value;//checked
};

double get_frequency(const wavelength& lambda_) 
{
    return 2*pi*SpeedOfLight/lambda_.value;//checked
};

double get_wavelength(const period& period_) 
{
    return SpeedOfLight*period_.value;//checked
};


//using previous three equatioons
double get_period(const wavelength& lambda_) 
{
    return lambda_.value/SpeedOfLight;
};

double get_wavelength(const frequency& omega_) 
{
    return 2*pi*SpeedOfLight/(omega_.value);
};

double get_frequency(const period& period_) 
{
    return 2.*pi/period_.value;
};




class Wave{
    //properties for something that is described by E(t) = Sin(omega*(t-t0)+phi)
    private:
        double Period = 0;
        double Omega = 0;
        double Lambda = 0;
        double InitialTime = 0;
        double Phase = 0;

    public:
        Wave(){};
        double operator()(const double& time);
        void set_Period(const double& Period_);
        void set_Omega(const double& Omega_);
        void set_Lambda(const double& WaveLength_);
        void set_InitialTime(const double& InitialTime_);
        void set_Phase(const double& Phase_);

        double get_Period();
        double get_Omega();
        double get_Lambda();
        double get_InitialTime();
        double get_Phase();
};


class Envelope{
    //properties for something that is described by E(t) = Sin^2((t-t0)/T) for 0<t-t0<T
    private:
        double Duration = 0;
        double InitialTime = 0;
    
    public:
        double operator()(const double& time);
        void set_Duration(const double& Duration_);
        void set_InitialTime(const double& InitialTime_);
};


class Laser{
    private:    
        std::vector<Coordinate<k>> Values;
        std::shared_ptr<TimeGrid> timegrid;

        Wave PlaneWave;
        Envelope envelope;

        double NumberOfCycles = 0;
        double Intensity = 0;
        double Amplitude = 0;
        Coordinate<k> Polarization;//you need this only in cartesian coordinates.

    public:
        Laser(){};
        auto operator()(const double& Time)
                {return (Amplitude*envelope(Time)*PlaneWave(Time))*Polarization;};
        //setters

        //only one of the following three calculates the others
        void set_Period(const double& Period_, const Unit& InputUnit) 
        {
            auto Period = Convert(Period_, InputUnit, AuTime); 
            PlaneWave.set_Period(Period);
        };
        void set_Omega(const double& Omega_)//, const Unit& InputUnit) 
        {
            //auto Omega = Convert(Omega_, InputUnit, AuTime); 
            PlaneWave.set_Omega(Omega_);
        };
        void set_Lambda(const double& WaveLength_, const Unit& InputUnit) 
        {
            auto lambda = Convert(WaveLength_, InputUnit, AuLength);
            PlaneWave.set_Lambda(lambda);
        };
        ////////////////////////////////////////////////////////
        void set_InitialTime(const double& InitialTime_, const Unit& InputUnit) 
        {
            PlaneWave.set_InitialTime(InitialTime_);  
            envelope.set_InitialTime(InitialTime_); 
        };
        void set_Phase(const double& Phase_)
        {
            PlaneWave.set_Phase(Phase_);
        };
        void set_NumberOfCycles(const double& NumberOfCycles_)
        {
            NumberOfCycles = NumberOfCycles_; 
            auto Duration = NumberOfCycles*PlaneWave.get_Period();
            std::cout << "Number of Cycles "<< NumberOfCycles << " Duration " << Convert(Duration, AuTime, FemtoSeconds) << std::endl;
            envelope.set_Duration(Duration);
        }
        void set_Intensity(const double& Intensity_, const Unit& InputUnit)
        {
            auto Intensity = Convert(Intensity_, InputUnit, AuIntensity);
            Amplitude = sqrt(Intensity);
        }


        void set_Polarization(const Coordinate<k>& Polarization_)
        {
            Polarization = Polarization_/Polarization_.norm(); 
        };
        void set_TimeGrid(const TimeGrid& TimeGrid_)
                {timegrid = std::make_shared<TimeGrid>(TimeGrid_);};

        void print_info()
        {
            std::cout << "********** Laser info **********\n";
            std::cout << "* frequency:  ";
            std::cout << std::setw(10) << std::setprecision(4) << Convert(PlaneWave.get_Omega(), AuEnergy, ElectronVolt) << " eV" << std::endl;
            std::cout << "* wavelength: ";
            std::cout << std::setw(10) << std::setprecision(4) << Convert(PlaneWave.get_Lambda(), AuLength, NanoMeters) << " nm" << std::endl;
            std::cout << "* Period:     ";
            std::cout << std::setw(10) << std::setprecision(4) << Convert(PlaneWave.get_Period(), AuTime, FemtoSeconds) << " fs" << std::endl;
            std::cout << "* Intensity:  ";
            std::cout << std::setw(10) << std::setprecision(4) << Amplitude*Amplitude << " au" << std::endl;
        }

};


class SetOfLaser
{
    private:
        std::vector<Laser> LaserArray;
    public:
        SetOfLaser(){};
        Laser& operator[](const int& i)
        {
            return const_cast<Laser&>(static_cast<const SetOfLaser&>(*this)[i]);
        };
        const Laser& operator[](const int& i) const
        {
            return LaserArray[i];
        };

        auto operator()(const double& Time)
        {
            Coordinate<k> result({0,0,0});
            for(auto& laser : LaserArray){
                result += laser(Time);
            }
            return result;
        };
};
#include "Laser_definitions.hpp"

#endif