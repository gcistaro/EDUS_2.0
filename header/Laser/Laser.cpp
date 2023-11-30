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

    if(t < 0 || t > Duration){
        return 0.;
    }
    return std::pow(std::sin(t/Duration), 2.);
}

void Envelope::set_Duration(const double& Duration_)
{
    Duration = Duration_;
}

void Envelope::set_InitialTime(const double& InitialTime_)
{
    InitialTime = InitialTime_;
}
