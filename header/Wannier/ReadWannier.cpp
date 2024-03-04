#include "StreamFile.hpp"
//each function reads part of wannier "_tb.dat" file

int ParseWannier_Degeneracies(const auto& LineIterator_Begin, auto& Degeneracy)
{  
    auto LineIterator_End = LineIterator_Begin;
    while((*LineIterator_End).size() > 0){
        for(auto& deg : (*LineIterator_End)){
            Degeneracy.push_back(std::atof(deg.c_str()));
        }
        ++LineIterator_End;
    }
    return (LineIterator_End - LineIterator_Begin);
}


int ParseWannier_MatrixElement(const auto& LineIterator, auto* Matrix_, const auto& NumberOfBands)
{
    assert( (*LineIterator).size() == 4 );
    auto m = std::atoi((*LineIterator)[0].c_str())-1;
    auto n = std::atoi((*LineIterator)[1].c_str())-1;
    assert( m<NumberOfBands && n<NumberOfBands );
    
    *(Matrix_+n+NumberOfBands*m) = std::atof((*LineIterator)[2].c_str()) + im*std::atof((*LineIterator)[3].c_str());
    return 1;
}

int ParseWannier_MatrixElement(const auto& LineIterator, auto* Matrix0,
auto* Matrix1, auto* Matrix2, const auto& NumberOfBands)
{
    assert( (*LineIterator).size() == 8 );

    auto m = std::atoi((*LineIterator)[0].c_str())-1;
    auto n = std::atoi((*LineIterator)[1].c_str())-1;
    assert( m<NumberOfBands && n<NumberOfBands );
    
    *(Matrix0+n+NumberOfBands*m) = std::atof((*LineIterator)[2].c_str()) + im*std::atof((*LineIterator)[3].c_str());
    *(Matrix1+n+NumberOfBands*m) = std::atof((*LineIterator)[4].c_str()) + im*std::atof((*LineIterator)[5].c_str());
    *(Matrix2+n+NumberOfBands*m) = std::atof((*LineIterator)[6].c_str()) + im*std::atof((*LineIterator)[7].c_str());
    return 1;
}

int ParseWannier_Matrix(const auto LineIterator_begin, auto* R, auto* Matrix_, const auto& NumberOfBands)
{
    auto LineIterator_aux = LineIterator_begin;
    *(R)   = std::atof((*LineIterator_aux)[0].c_str());
    *(R+1) = std::atof((*LineIterator_aux)[1].c_str());
    *(R+2) = std::atof((*LineIterator_aux)[2].c_str());

    LineIterator_aux++;
    while((*LineIterator_aux).size()>0){
        int SizeOfMatrixElement = ParseWannier_MatrixElement(LineIterator_aux, Matrix_, NumberOfBands);
        LineIterator_aux+=SizeOfMatrixElement;
    }
    return LineIterator_aux-LineIterator_begin;
}


int ParseWannier_Matrix(const auto LineIterator_begin, auto* R, auto* Matrix0, auto* Matrix1, auto* Matrix2, const auto& NumberOfBands)
{
    auto LineIterator_aux = LineIterator_begin;
    *(R)   = std::atof((*LineIterator_aux)[0].c_str());
    *(R+1) = std::atof((*LineIterator_aux)[1].c_str());
    *(R+2) = std::atof((*LineIterator_aux)[2].c_str());
    
    LineIterator_aux++;
    while((*LineIterator_aux).size()>0){
        int SizeOfMatrixElement = ParseWannier_MatrixElement(LineIterator_aux, Matrix0, Matrix1, Matrix2, NumberOfBands);
        LineIterator_aux+=SizeOfMatrixElement;
    }
    return LineIterator_aux-LineIterator_begin;
}

int ParseWannier_Hamiltonian(const auto LineIterator_begin, auto& Rmesh, auto& H, const auto& NumberOfBands)
{
    auto LineIterator_aux = LineIterator_begin;
    int index = 0;
    while((*(LineIterator_aux+1)).size() == 4){//the things in while goes to another block and '+1' skips the Rvector
        int SizeOfMatrix = ParseWannier_Matrix(LineIterator_aux, &Rmesh(index,0), &H(index,0,0), NumberOfBands);
        LineIterator_aux += SizeOfMatrix+1;
        ++index;
    }
    return LineIterator_aux-LineIterator_begin;
}

int ParseWannier_PositionOperator(const auto& LineIterator_begin, auto& Rmesh, auto& r, const auto& NumberOfBands)
{
    auto LineIterator_aux = LineIterator_begin;
    int index = 0;
    while((*(LineIterator_aux+1)).size() == 8){//the things in while goes to another block and '+1' skips the Rvector
        int SizeOfMatrix = ParseWannier_Matrix(LineIterator_aux, &Rmesh(index,0), &(r[0](index,0,0)),
                                               &(r[1](index,0,0)), &(r[2](index,0,0)), NumberOfBands);
        LineIterator_aux += SizeOfMatrix+1;
        ++index;
    }
    return LineIterator_aux-LineIterator_begin;
}

void ParseWannier(const std::string& FileNameTB, auto& NumberOfBands, auto& NumberOfRpoints,
                  auto& UnitCell, auto& Degeneracy, auto& Rmesh, auto& H, auto& r)
{
    auto file_content = ReadFile(FileNameTB);
    
    //read unit cell
    UnitCell.initialize({3,3});
    auto LineIterator = file_content.begin()+1;
    assert((*LineIterator).size() == 3);
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            UnitCell(i,j) = std::atof((*LineIterator)[j].c_str());
        }
        ++LineIterator;
    }
    
    NumberOfBands   = std::atoi((*(file_content.begin() + 4))[0].c_str());
    NumberOfRpoints = std::atoi((*(file_content.begin() + 5))[0].c_str());

    LineIterator = file_content.begin()+6;
    auto Line_EndDegeneracy = ParseWannier_Degeneracies(LineIterator, Degeneracy);    
    assert(Degeneracy.size() == NumberOfRpoints);

    ////now we are at the blank line before the beginning of matrix elements
    LineIterator+=Line_EndDegeneracy;
    while((*LineIterator).size() == 0){
        ++LineIterator;
    }

    //Read Hamiltonian Blocks
    H.initialize({NumberOfRpoints, NumberOfBands, NumberOfBands});
    Rmesh.initialize({NumberOfRpoints, 3});
    r[0].initialize({NumberOfRpoints, NumberOfBands, NumberOfBands});
    r[1].initialize({NumberOfRpoints, NumberOfBands, NumberOfBands});
    r[2].initialize({NumberOfRpoints, NumberOfBands, NumberOfBands});
    
    auto Line_endHamiltonian = ParseWannier_Hamiltonian(LineIterator, Rmesh, H, NumberOfBands);
    LineIterator += Line_endHamiltonian;


    auto Line_endPositionOperator = ParseWannier_PositionOperator(LineIterator, Rmesh, r, NumberOfBands);
    
}


