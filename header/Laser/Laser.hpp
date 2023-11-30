#ifndef LASER_HPP
#define LASER_HPP

//#include <math.h>

#include "../Constants.hpp"
#include "../ConvertUnits.hpp"
#include "../Geometry/Geometry.hpp"
#include "../TimeGrid/TimeGrid.hpp"

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
    return 2.*pi/omega_.value;
};

double get_frequency(const wavelength& lambda_) 
{
    return SpeedOfLight/lambda_.value;
};

double get_wavelength(const period& period_) 
{
    return SpeedOfLight*period_.value/(2*pi);
};


//using previous three equatioons
double get_period(const wavelength& lambda_) 
{
    return 2.*pi*lambda_.value/SpeedOfLight;
};

double get_wavelength(const frequency& omega_) 
{
    return SpeedOfLight/omega_.value;
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
        std::vector<Vector<k>> Values;
        std::shared_ptr<TimeGrid> timegrid;

        Wave PlaneWave;
        Envelope envelope;

        double NumberOfCycles = 0;
        double Intensity = 0;
        double Amplitude = 0;
        Vector<k> Polarization;//you need this only in cartesian coordinates.

    public:
        Laser(){};
        auto operator()(const double& Time)
                {return (Amplitude*envelope(Time)*PlaneWave(Time))*Polarization;};
        //setters

        //only one of the following three calculates the others
        void set_Period(const double& Period_) 
                {PlaneWave.set_Period(Period_);};
        void set_Omega(const double& Omega_)
                {PlaneWave.set_Omega(Omega_);};
        void set_Lambda(const double& WaveLength_)
                {PlaneWave.set_Lambda(WaveLength_);};
        ////////////////////////////////////////////////////////
        void set_InitialTime(const double& InitialTime_)
                {PlaneWave.set_InitialTime(InitialTime_);  envelope.set_InitialTime(InitialTime_); };
        void set_Phase(const double& Phase_)
                {PlaneWave.set_Phase(Phase_);};
        void set_NumberOfCycles(const double& NumberOfCycles_)
                {NumberOfCycles = NumberOfCycles_; auto Duration = NumberOfCycles*PlaneWave.get_Period();
                 envelope.set_Duration(Duration);}
        void set_Amplitude(const double& Amplitude_)
                {Amplitude = Amplitude_;}


        void set_Polarization(const Vector<k>& Polarization_)
                {Polarization = Polarization_; };
        void set_TimeGrid(const TimeGrid& TimeGrid_)
                {timegrid = std::make_shared<TimeGrid>(TimeGrid_);};
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
            Vector<k> result({0,0,0});
            for(auto& laser : LaserArray){
                result += laser(Time);
            }
            return result;
        };
};
#include "Laser.cpp"

#endif