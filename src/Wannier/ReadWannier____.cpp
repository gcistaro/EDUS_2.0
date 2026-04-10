#include "StreamFile.hpp"
#include "ReadWannier.hpp"
//each function reads part of wannier "_tb.dat" file

int ParseWannier_Degeneracies(const std::vector<std::vector<std::string>>::iterator& LineIterator_Begin, 
                              std::vector<int>& Degeneracy)
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


int ParseWannier_MatrixElement(const std::vector<std::vector<std::string>>::iterator& LineIterator, 
                                std::complex<double>* Matrix_, const int& NumberOfBands)
{
    assert( (*LineIterator).size() == 4 );
    auto m = std::atoi((*LineIterator)[0].c_str())-1;
    auto n = std::atoi((*LineIterator)[1].c_str())-1;
    assert( m<NumberOfBands && n<NumberOfBands );
    
    *(Matrix_+n+NumberOfBands*m) = std::atof((*LineIterator)[2].c_str()) + im*std::atof((*LineIterator)[3].c_str());
    return 1;
}

int ParseWannier_MatrixElement(const std::vector<std::vector<std::string>>::iterator& LineIterator, 
                                std::complex<double>* Matrix0, std::complex<double>* Matrix1, 
                                std::complex<double>* Matrix2, const int& NumberOfBands)
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

int ParseWannier_Matrix(const std::vector<std::vector<std::string>>::iterator LineIterator_begin, 
                        double* R, std::complex<double>* Matrix_, const int& NumberOfBands)
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


int ParseWannier_Matrix(const std::vector<std::vector<std::string>>::iterator LineIterator_begin, 
                        const std::vector<std::vector<std::string>>::iterator& LineIterator_end,
                        double* R, std::complex<double>* Matrix0, std::complex<double>* Matrix1, 
                        std::complex<double>* Matrix2, const int& NumberOfBands)
{
    auto LineIterator_aux = LineIterator_begin;
    *(R)   = std::atof((*LineIterator_aux)[0].c_str());
    *(R+1) = std::atof((*LineIterator_aux)[1].c_str());
    *(R+2) = std::atof((*LineIterator_aux)[2].c_str());
    
    LineIterator_aux++;
    while(( (LineIterator_aux != LineIterator_end) &&
            (*(LineIterator_aux)).size() == 8)  ) {//the things in while goes to another block and '+1' skips the Rvector
        
        int SizeOfMatrixElement = ParseWannier_MatrixElement(LineIterator_aux, Matrix0, Matrix1, Matrix2, NumberOfBands);
        LineIterator_aux+=SizeOfMatrixElement;
    }
    return LineIterator_aux-LineIterator_begin;
}

int ParseWannier_Hamiltonian(const std::vector<std::vector<std::string>>::iterator LineIterator_begin, 
                            mdarray<double,2>& Rmesh, mdarray<std::complex<double>, 3>& H, 
                            const int& NumberOfBands)
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

int ParseWannier_PositionOperator(const std::vector<std::vector<std::string>>::iterator& LineIterator_begin,
                 const std::vector<std::vector<std::string>>::iterator& LineIterator_end,
                 mdarray<double,2>& Rmesh, std::array<mdarray<std::complex<double>,3>, 3>& r, 
                 const int& NumberOfBands)
{
    auto LineIterator_aux = LineIterator_begin;
    int index = 0;
    while(( (LineIterator_aux+1 < LineIterator_end) && 
            (*(LineIterator_aux+1)).size() == 8)  ) {//the things in while goes to another block and '+1' skips the Rvector
        int SizeOfMatrix = ParseWannier_Matrix(LineIterator_aux, LineIterator_end, &Rmesh(index,0), &(r[0](index,0,0)),
                                               &(r[1](index,0,0)), &(r[2](index,0,0)), NumberOfBands);
        LineIterator_aux += SizeOfMatrix+1;
        ++index;
    }
    return LineIterator_aux-LineIterator_begin;
}

void ParseWannier(const std::string& FileNameTB, int& NumberOfBands, int& NumberOfRpoints,
                  mdarray<double,2>& UnitCell, std::vector<int>& Degeneracy, mdarray<double,2>& Rmesh, 
                  mdarray<std::complex<double>, 3>& H, std::array<mdarray<std::complex<double>,3>, 3>& r)
{
    auto file_content = ReadFile(FileNameTB);
    
    int iline = 1; //index to iterate over the lines 

    /* read unit cell */
    UnitCell.initialize({3,3});
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            UnitCell(i,j) = std::atof(file_content[iline+i][j].c_str());
        }
    }
    
    /* read number of bands and number of R points */
    NumberOfBands   = std::atoi(file_content[4][0].c_str());
    NumberOfRpoints = std::atoi(file_content[5][0].c_str());

    /* read degeneracies */
    Degeneracy.resize(NumberOfRpoints);
    double idegeneracy = 0;
    iline = 6;
    int iword = 0; //counter over numbers in a line
    do {
        if ( iword == file_content[iline].size() ) {
            iword = 0;
            iline++;
        }
        Degeneracy[idegeneracy] = std::atoi(file_content[iline][iword].c_str());
        idegeneracy++;
        iword++;
    }
    while( idegeneracy < Degeneracy.size());

    iline++;
    if( file_content[iline].size() > 0 ) {
        std::runtime_error("No blank line in tb file after degeneracies.\n");
    }

    /* read Hamiltonian */
    int iR = 0;
    Rmesh.initialize({NumberOfRpoints,3});
    H.initialize({NumberOfRpoints, NumberOfBands, NumberOfBands});

    while ( iR < NumberOfRpoints ) {
        iline++;
        if( file_content[iline].size() != 3 ) {
            std::runtime_error("Expected R vector at line but not found.\n");
        }

        /* read R vector */
        for( auto& ix : {0,1,2} ) {
            Rmesh(iR, ix) = std::atoi(file_content[iline][ix].c_str());
        }
        /* read matrix at R */
        iline++;
        while ( file_content[iline].size() > 0 ) {
            auto irow_fort = std::atoi(file_content[iline][0].c_str());
            auto icol_fort = std::atoi(file_content[iline][1].c_str());
            auto irow_c = irow_fort - 1;
            auto icol_c = icol_fort - 1;
            H(iR, irow_c, icol_c) = std::atof(file_content[iline][2].c_str()) +
                                 im*std::atof(file_content[iline][3].c_str());
            iline++;
        }
        iR++;

    }
std::cout << "hamiltonain read " << std::endl;
    /* read r operator */
    r[0].initialize({NumberOfRpoints, NumberOfBands, NumberOfBands});
    r[1].initialize({NumberOfRpoints, NumberOfBands, NumberOfBands});
    r[2].initialize({NumberOfRpoints, NumberOfBands, NumberOfBands});

    iR = 0;
    while ( iR < NumberOfRpoints ) {
        std::cout << "ir" << iR << std::endl;
        iline++;
        if( file_content[iline].size() != 3 ) {
            std::runtime_error("Expected R vector at line but not found.\n");
        }
        /* read R vector */
        for( auto& ix : {0,1,2} ) {
            Rmesh(iR, ix) = std::atoi(file_content[iline][ix].c_str());
        }
        /* read matrix at R */
        iline++;
        while ( file_content.size() < iline && file_content[iline].size() > 0 ) {
            std::cout << file_content[iline][0] << std::endl;
            auto irow_fort = std::atoi(file_content[iline][0].c_str());
            auto icol_fort = std::atoi(file_content[iline][1].c_str());
            auto irow_c = irow_fort - 1;
            auto icol_c = icol_fort - 1;
            r[0](iR, irow_c, icol_c) = std::atof(file_content[iline][2].c_str()) +
                             im*std::atof(file_content[iline][3].c_str());
            r[1](iR, irow_c, icol_c) = std::atof(file_content[iline][4].c_str()) +
                             im*std::atof(file_content[iline][5].c_str());
            r[2](iR, irow_c, icol_c) = std::atof(file_content[iline][6].c_str()) +
                             im*std::atof(file_content[iline][7].c_str());
            iline++;
        }
        iR++;
    }
std::cout << "r read " << std::endl;

}


