//each function reads part of wannier "_tb.dat" file
#ifndef READWANNIER_HPP
#define READWANNIER_HPP


#include <stdlib.h>//for atof and atoi

#include "../Constants.hpp"
#include "../mdContainers/mdContainers.hpp"
#include "../StreamFile.hpp"

//each function reads part of wannier "_tb.dat" file

int ParseWannier_Degeneracies(const auto& LineIterator_Begin, auto& Degeneracy);

int ParseWannier_MatrixElement(const auto& LineIterator, auto* Matrix_, const auto& NumberOfBands);
int ParseWannier_MatrixElement(const auto& LineIterator, auto* Matrix0, auto* Matrix1, auto* Matrix2, const auto& NumberOfBands);

int ParseWannier_Matrix(const auto LineIterator_begin, auto* R, auto* Matrix_, const auto& NumberOfBands);
int ParseWannier_Matrix(const auto LineIterator_begin, auto* R, auto* Matrix0, auto* Matrix1, auto* Matrix2, const auto& NumberOfBands);

int ParseWannier_Hamiltonian(const auto LineIterator_begin, auto& Rmesh, auto& H, const auto& NumberOfBands);
int ParseWannier_PositionOperator(const auto& LineIterator_begin, auto& temporary_Rvectors, auto& r, const auto& NumberOfBands);

void ParseWannier(const std::string& FileNameTB, auto& NumberOfBands, auto& NumberOfRpoints,
                  auto& UnitCell, auto& Degeneracy, auto& Rmesh, auto& H, auto& r);



#include "ReadWannier.cpp"

#endif
