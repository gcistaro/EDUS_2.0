#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>

#include "../Constants.hpp"


class TimeGrid{
    private:
        double ResolutionTime=0.;
        size_t NumberOfSteps=0;
        double InitialTime=0.;
        double FinalTime=0.;
        std::vector<double> timegrid;

    public:
        TimeGrid(){};
        TimeGrid(const double& InitialTime_, const double& FinalTime_, const double& ResolutionTime_)
        : ResolutionTime(ResolutionTime_), InitialTime(InitialTime_), FinalTime(FinalTime_)
        {
            this->initialize();
        }        

        void initialize()
        {
            assert(FinalTime-InitialTime > ResolutionTime && ResolutionTime > threshold);
            
            NumberOfSteps = size_t(FinalTime-InitialTime)/ResolutionTime;
            timegrid.resize(NumberOfSteps);

            for(size_t step=0; step<NumberOfSteps; step++){
                timegrid[step] = InitialTime + step*ResolutionTime;
            }
        }

        friend std::ostream& operator<<(std::ostream& os, TimeGrid& tg)
        {
            for(auto& time : tg.timegrid){
                os << std::setw(12) << std::setprecision(5) << time << std::endl; 
            }
            return os;
        }

};