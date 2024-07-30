#include "PrintWannier.hpp"

void PrintWannier(const std::string& FileName, const int& NumberOfBands, const int& NumberOfRpoints, 
                   const mdarray<double,2>& UnitCell, const std::vector<int>& Degeneracy, const mdarray<double,2>& Rmesh, 
                   const mdarray<std::complex<double>,3>& H, const std::array<mdarray<std::complex<double>,3>, 3>& r)
{
    std::ofstream OutputFile;
    OutputFile.open(FileName.c_str());

    OutputFile << "#Created with NEGF program for checking purposes" << std::endl;
    
    //unitcell
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            OutputFile << std::setw(10) << UnitCell(i,j);
        }
        OutputFile << std::endl;
    }
    //end unitcell


    OutputFile << NumberOfBands;
    OutputFile << std::endl;
    OutputFile << NumberOfRpoints;

    //degeneracy
    int counter_deg = 0;
    for(auto& deg_ : Degeneracy){
        if(counter_deg % 15 ==0 ){
            OutputFile << std::endl;
        }
        OutputFile << std::setw(5) << deg_;
        ++counter_deg;
    }
    OutputFile << std::endl;
    OutputFile << std::endl;
    //end degeneracy

    //Hamiltonian
    for(int iR=0; iR<NumberOfRpoints; iR++){
        OutputFile << std::setw(5) << Rmesh(iR, 0);
        OutputFile << std::setw(5) << Rmesh(iR, 1);
        OutputFile << std::setw(5) << Rmesh(iR, 2);
        OutputFile << std::endl;
        for(int column=0; column<NumberOfBands; column++){
            for(int row=0; row<NumberOfBands; row++){
                OutputFile << std::setw(5) << row+1;    
                OutputFile << std::setw(5) << column+1;
                OutputFile << std::setw(10) << H(iR,row,column).real();
                OutputFile << std::setw(10) << H(iR,row,column).imag();
                OutputFile << std::endl;
            }
        }
        OutputFile << std::endl;
    }

    //Position operator
    for(int iR=0; iR<NumberOfRpoints; iR++){
        OutputFile << std::setw(5) << Rmesh(iR,0);
        OutputFile << std::setw(5) << Rmesh(iR,1);
        OutputFile << std::setw(5) << Rmesh(iR,2);
        OutputFile << std::endl;
        for(int column=0; column<NumberOfBands; column++){
            for(int row=0; row<NumberOfBands; row++){
                OutputFile << std::setw(5) << row+1;    
                OutputFile << std::setw(5) << column+1;
                OutputFile << std::setw(10) << r[0](iR,row,column).real();
                OutputFile << std::setw(10) << r[0](iR,row,column).imag();
                OutputFile << std::setw(10) << r[1](iR,row,column).real();
                OutputFile << std::setw(10) << r[1](iR,row,column).imag();
                OutputFile << std::setw(10) << r[2](iR,row,column).real();
                OutputFile << std::setw(10) << r[2](iR,row,column).imag();
                OutputFile << std::endl;
            }
        }
        OutputFile << std::endl;
    }

    //end Position operator

    OutputFile.close();


}

