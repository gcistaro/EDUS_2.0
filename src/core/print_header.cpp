#include <iomanip>
#include "print_header.hpp"
#include "GlobalFunctions.hpp"

void print_header()
{
    time_t now = time(0);
   
    // convert now to string form
    char* dt = ctime(&now);
    
    output::stars();
    output::print(std::string(21, ' '), "        ________   ________     ___    ___     ________");                   
    output::print(std::string(21, ' '), "       |  _____|| |  __   \\\\   |  ||  |  ||   /  _____||");                
    output::print(std::string(21, ' '), " ___   |  ||___   | ||  \\  \\\\  |  ||  |  ||   | ||____      ___");        
    output::print(std::string(21, ' '), "|___|  |   ___||  | ||   | ||  |  ||  |  ||   \\_____  \\\\   |___|");       
    output::print(std::string(21, ' '), "       |  ||____  | ||__/  //  |  \\\\_/   ||    _____\\  ||");                 
    output::print(std::string(21, ' '), "       |_______|| |_______//    \\_______//    |_______//");        
    output::print(" ");            
    output::stars();


    output::title("GIT RECAP");
    std::stringstream start;
    start << "Execution started: " << dt;
    output::print(start.str());

    output::print("Git hash:           *    ", git_hash);
    output::print("Git branch:         *    ", git_branchname);
    output::stars();
}
