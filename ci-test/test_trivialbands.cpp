#include "Simulation/Simulation.hpp"
#include "core/print_timing.hpp"
//OUTDATED!!! REMAKE IT!!! it was useful
/*
This test is constructed with simple Hamiltonian and simple eigenstates .
We use:

Energy[k](0) = -3/2 + Cos(2pi i k)
Energy[k](1) = +3/2 + Cos(2pi i k)

And eigenstates:
       | Exp[+2pi i k]   -Exp[+2pi i k] |
U[k] = |                                |
       | Exp[-2pi i k]    Exp[-2pi i k] |

So that 

         | Delta[0,R]       Delta[-2,R] |
Rho[R] = |                              |
         | Delta[+2,R]      Delta[0,R]  |
 
*/      
int main()
{            
    Simulation simulation("/home/gcistaro/EDUS/tb_models/2B_CosBand", 20.);
    //TODO: assert dm has the expected form
    print_timing(1);
}
