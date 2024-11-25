#include "Laser/Laser.hpp"

double get_period(const frequency& omega_) 
{
    return 2.*pi/omega_.value;//checked
}

double get_frequency(const wavelength& lambda_) 
{
    return 2*pi*SpeedOfLight/lambda_.value;//checked
}

double get_wavelength(const period& period_) 
{
    return SpeedOfLight*period_.value;//checked
}

//using previous three equatioons
double get_period(const wavelength& lambda_) 
{
    return lambda_.value/SpeedOfLight;
}

double get_wavelength(const frequency& omega_) 
{
    return 2*pi*SpeedOfLight/(omega_.value);
}

double get_frequency(const period& period_) 
{
    return 2.*pi/period_.value;
}


double Wave::operator()(const double& Time)
{
    return std::sin(Omega*(Time - InitialTime) + Phase);
}

void Wave::set_Period(const double& Period_)
{
    Period = Period_;
    Omega = get_frequency(period(Period));
    Lambda = get_wavelength(period(Period));
}

void Wave::set_Omega(const double& Omega_)
{
    Omega = Omega_;
    Period = get_period(frequency(Omega));
    Lambda = get_wavelength(frequency(Omega));
}

void Wave::set_Lambda(const double& WaveLength_)
{
    Lambda = WaveLength_;
    Period = get_period(wavelength(Lambda));
    Omega = get_frequency(wavelength(Lambda));
}

void Wave::set_InitialTime(const double& InitialTime_)
{
    InitialTime = InitialTime_;
}

void Wave::set_Phase(const double& Phase_)
{
    Phase = Phase_;
}

double Wave::get_Period()
{
    return this->Period;
}

double Wave::get_Omega()
{
    return this->Omega;
}

double Wave::get_Lambda()
{
    return this->Lambda;
}

double Wave::get_InitialTime()
{
    return this->InitialTime;
}
double Wave::get_Phase()
{
    return this->Phase;
}


double Envelope::operator()(const double& Time)
{
    double t = Time - InitialTime;

    if( std::abs(t) <= 1.e-07 || t > Duration){
        return 0.;
    }
    return std::pow(std::sin(pi*t/Duration), 2.);
}

void Envelope::set_Duration(const double& Duration_)
{
    Duration = Duration_;
}

double Envelope::get_Duration()
{
    return Duration;
}

double Envelope::get_InitialTime()
{
    return InitialTime;
}

double Envelope::get_FinalTime()
{
    return Duration + InitialTime;
}

void Envelope::set_InitialTime(const double& InitialTime_)
{
    InitialTime = InitialTime_;
}

Coordinate Laser::operator()(const double& Time)
{
    return (Amplitude*envelope(Time)*PlaneWave(Time))*Polarization;
}

Coordinate Laser::VectorPotential(const double& Time)
{
    //integral of the electric field with a - sign 
    if(Time < get_InitialTime() || Time > get_Duration()) {
        return Coordinate(0,0,0);
    }
    auto a = pi/envelope.Duration;
    auto b = PlaneWave.Omega;
    auto one_ = 2.*a-b;
    auto two_ = 2*a+b;
    auto three_ = b;
    auto& t0 = envelope.InitialTime;
    return -(.25*Amplitude*(-std::cos(one_*Time)/(one_)+std::cos(two_*Time)/(two_)-2.*std::cos(three_*Time)/(three_))
            -.25*Amplitude*(-std::cos(one_*t0)/(one_)+std::cos(two_*t0)/(two_)-2.*std::cos(three_*t0)/(three_)))*Polarization;
}

//only one of the following three calculates the others
void Laser::set_Period(const double& Period_, const Unit& InputUnit) 
{
    auto Period = Convert(Period_, InputUnit, AuTime); 
    PlaneWave.set_Period(Period);
}

void Laser::set_Omega(const double& Omega_, const Unit& InputUnit)//, const Unit& InputUnit) 
{
    auto Omega = Convert(Omega_, InputUnit, AuEnergy);
    PlaneWave.set_Omega(Omega);
}

void Laser::set_Lambda(const double& WaveLength_, const Unit& InputUnit) 
{
    auto lambda = Convert(WaveLength_, InputUnit, AuLength);
    PlaneWave.set_Lambda(lambda);
}
////////////////////////////////////////////////////////

void Laser::set_InitialTime(const double& InitialTime_, const Unit& InputUnit) 
{
    PlaneWave.set_InitialTime(InitialTime_);  
    envelope.set_InitialTime(InitialTime_); 
}

void Laser::set_Phase(const double& Phase_)
{
    PlaneWave.set_Phase(Phase_);
}
void Laser::set_NumberOfCycles(const double& NumberOfCycles_)
{
    NumberOfCycles = NumberOfCycles_; 
    auto Duration = NumberOfCycles*PlaneWave.get_Period();
    envelope.set_Duration(Duration);
}

void Laser::set_Intensity(const double& Intensity_, const Unit& InputUnit)
{
    auto Intensity = Convert(Intensity_, InputUnit, AuIntensity);
    Amplitude = sqrt(Intensity);
}

void Laser::set_Polarization(const Coordinate& Polarization_)
{
    Polarization = Polarization_/Polarization_.norm(); 
}

void Laser::set_TimeGrid(const TimeGrid& TimeGrid_)
{
    timegrid = std::make_shared<TimeGrid>(TimeGrid_);
}

double Laser::get_Duration()
{
    return envelope.get_Duration();
}

double Laser::get_InitialTime()
{
    return envelope.get_InitialTime();
}

double Laser::get_FinalTime()
{
    return envelope.get_FinalTime();
}

double Laser::get_NumberOfCycles()
{
    return NumberOfCycles;
}

double Laser::get_Omega() 
{
    return PlaneWave.get_Omega();
}

double Laser::get_Lambda() 
{
    return PlaneWave.get_Lambda();
}

void Laser::print_info()
{
    std::cout << "****************************************************       LASER INFO        *************************************************\n";
    std::cout << "*   Frequency:          *     ";
    std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << PlaneWave.get_Omega();
    std::cout << std::left << std::setw(15) << " a.u.";
    std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << Convert(PlaneWave.get_Omega(), AuEnergy, ElectronVolt);
    std::cout << std::left << std::setw(60) << " eV" << "*\n";
    std::cout << "*   Wavelength:         *     ";
    std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << PlaneWave.get_Lambda();
    std::cout << std::left << std::setw(15) << " a.u.";
    std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << Convert(PlaneWave.get_Lambda(), AuLength, NanoMeters);
    std::cout << std::left << std::setw(60) << " nm" << "*\n";
    std::cout << "*   Period:             *     ";
    std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << PlaneWave.get_Period();
    std::cout << std::left << std::setw(15) << " a.u.";
    std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << Convert(PlaneWave.get_Period(), AuTime, FemtoSeconds);
    std::cout << std::left << std::setw(60) << " fs" << "*\n";
    std::cout << "*   Intensity:          *     ";
    std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << Amplitude*Amplitude;
    std::cout << std::left << std::setw(15) << " a.u.";
    std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << Convert(Amplitude*Amplitude, AuIntensity, Wcm2);
    std::cout << std::left << std::setw(60) << " Wcm2" << "*\n";
    std::cout << "******************************************************************************************************************************\n";
}

