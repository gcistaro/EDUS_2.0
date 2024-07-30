#ifndef PRINTHEADER_HPP
#define PRINTHEADER_HPP

#include <ctime>
#include <iostream>

#define watch(x)  (#x)
#include "core/githash.hpp"

void print_header()
{
    time_t now = time(0);
   
    // convert now to string form
    char* dt = ctime(&now);
    
    std::cout << "Execution started: " << dt << std::endl;
    std::cout << "*********************************************************************************\n";
    std::cout << "*      _________           ________   ________     ___    ___     ________      *\n";
    std::cout << "*     |    _   ||         |  _____|| |  __   \\\\   |  ||  |  ||   /  _____||     *\n";
    std::cout << "*     |   |_|  ||   ___   |  ||___   | ||  \\  \\\\  |  ||  |  ||   | ||____       *\n";
    std::cout << "*     |        //  |___|  |   ___||  | ||   | ||  |  ||  |  ||   \\_____  \\\\     *\n";
    std::cout << "*     |   |\\   \\\\         |  ||____  | ||__/  //  |  \\\\_/   ||    _____\\  ||    *\n";
    std::cout << "*     |___| \\__||         |_______|| |_______//    \\_______//    |_______//     *\n";
    std::cout << "*                                                                               *\n";
    std::cout << "*********************************************************************************\n";
    std::cout << "Git hash:            " << git::hash << "\n";
    std::cout << "Git branch:          " << git::branchname << "\n";
    std::cout << "*********************************************************************************\n";
}

#endif