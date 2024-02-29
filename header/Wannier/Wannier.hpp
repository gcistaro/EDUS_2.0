#ifndef WANNIER_HPP
#define WANNIER_HPP

#include <string>
#include <iostream>
#include <sstream>

#include <vector>
#include "core/profiler.hpp"
#include "mdContainers/mdContainers.hpp"
#include "ReadWannier.hpp"
#include "PrintWannier.hpp"


class Wannier{
    private:
        std::vector<int> Degeneracy;
        mdarray<std::complex<double>, 3> H;
        std::array<mdarray<std::complex<double>,3>, 3> r;
        mdarray<double,2> Rmesh;
        mdarray<double,2> UnitCell;

        size_t NumberOfBands;
        size_t NumberOfRpoints;
    public: 
        Wannier(const std::string& FileName);
        void Print(const std::string& FileName);

        friend class Material;
};


Wannier::Wannier(const std::string& FileName)
{
    PROFILE("Wannier");
    std::stringstream TotalFileName;
    TotalFileName << FileName << "_tb.dat";    
    ParseWannier(TotalFileName.str(), NumberOfBands, NumberOfRpoints,
                  UnitCell, Degeneracy, Rmesh, H, r);
}

void Wannier::Print(const std::string& FileName)
{
    PrintWannier(FileName, NumberOfBands, NumberOfRpoints, 
                UnitCell, Degeneracy, Rmesh, H, r);
}


#endif



/*
ABOUT THE DEGENERACIES:
http://www.democritos.it/pipermail/wannier/2014-January/000757.html
If your ab-initio k-point grid contains num_kpts = N_1*N_2*N_3 k-points,
the Wannier functions are periodic over a supercell with dimensions
A_i=N_i*a_i, where {a_i} are the crystal-cell lattice vectors, and
i=1,2,3. The vectors {A_i} are the primitive translations of a
superlattice of vectors RR.

The optimal choice of R-vectors to be used in the Fourier sums are the
ones which lie within a Wigner-Seitz (W-S) supercell centered at the
origin: they are closer to RR=0 than to any other superlattice vectors RR.

You are right to expect that there should be exactly num_kpts such
R-vectors (i.e., nrpts = num_kpts). But there is a subtlety: sometimes you
will find R-vectors which are equidistant to RR=0 and to one or several
adjacent RR-vectors. When that happens, the R-vector in question is added
to your set, but the code also records, in the variable ndegen, the number
of W-S cells which "share" this R-vector. At the end you arrive at a set
of nrpts >= num_kpts R-vectors, such that the "sum rule"

sum_R 1/ndegen(R) = num_kpts

is satisfied, with the sum over the nrpts R-vectors. Of course, ndegen=1
for "interior" R-vectors (those which are not at W-S supercell boundary).
If all R-vectors are interior, then nrpts=num_kpts.
Suppose you have calculated the Hamiltonian matrix H_mn(R)=<0m|H|Rn> in
the Wannier basis. The Fourier transform H_mn(k) is obtained as

H_mn(k) = sum_R e^{+ik.R}*H_mn(R)/ndegen(R)

[note the weighting factor 1/ndegen(R)].
*/