#include "kGradient.hpp"

#include "omp.h"//TODO:: FIX ME PUT ME IN THE RIGHT PLACE!!


kGradient::kGradient(const MeshGrid& kmesh__) 
{
    this->initialize(kmesh__);
}

void kGradient::initialize(const MeshGrid& kmesh__)
{
    assert(kmesh__.get_type() == cube);
    kmesh = std::make_shared<MeshGrid>( get_GammaCentered_grid(kmesh__) );
    initialize();
}



void kGradient::initialize()
{
    //evaluates weights, number of shells, k+b indices.
    if(!kmesh) {
        return;
    }
    ikshell = SortInShells(*kmesh);
    Calculate_nshellsAndweights(nshells, Weight, *kmesh, ikshell);
    ikpb = Find_kpb(*kmesh, ikshell);
}

std::vector<std::vector<int>> SortInShells(const MeshGrid& kmesh)
{
    std::vector<std::vector<int>> ikshell; //[ishell][ik] contains all ik with norm ishell
    //find all norms of k vectors
    std::vector<double> knorms; //array with norm of k vectors (not repeated)
    for( auto& k_ : kmesh.get_mesh() ) {
        auto k_norm = k_.norm();
        auto it = std::find_if(knorms.begin(), knorms.end(), 
                               [&k_norm](const auto& norm) { return std::abs( norm - k_norm ) < threshold; });
        bool found = ( it != knorms.end() );
        if( !found ) {
            knorms.push_back( k_norm );
        }                 
    }
    std::sort( knorms.begin(), knorms.end() );
    ikshell.resize(knorms.size());
    //find indices of each shell
    for( size_t ik=0; ik< kmesh.get_TotalSize(); ++ik ) {
        //find what shell k_ belongs to
        auto k_norm = kmesh[ik].norm();
        auto it = std::find_if(knorms.begin(), knorms.end(), 
                               [&k_norm](const auto& norm) { return std::abs( norm - k_norm ) < 1.e-07; });
        assert(it != knorms.end());
        auto ishell = it - knorms.begin();
        //push back the index in ikshell
        ikshell[ishell].push_back(ik);
    }
    //remove 0 from kshells 
    assert( std::abs(knorms[0]) < threshold );
    ikshell.erase(ikshell.begin());
    knorms.erase(knorms.begin());

    //recap
    std::ofstream of("Cartesian.txt");
    for(int ishell=0; ishell<ikshell.size(); ishell++) {
        for(int ik=0; ik<ikshell[ishell].size(); ik++) {
            of << ishell << " " <<knorms[ishell] << " " <<  ikshell[ishell][ik] << " " << kmesh[ikshell[ishell][ik]].get("Cartesian");//get(LatticeVectors(k));
        }
    }
    of.close();

    return ikshell;
}


void Calculate_nshellsAndweights(int& nshells, mdarray<double,1>& Weight, 
                                 const MeshGrid& kmesh, const std::vector<std::vector<int>>& ikshell)
{    
    auto q = Matrix<double>({6,1});

    q(0,0) = ( kmesh.get_Size()[0] > 1 ? 1 : 0 ) ;     
    q(1,0) = ( kmesh.get_Size()[1] > 1 ? 1 : 0 ) ;   
    q(2,0) = ( kmesh.get_Size()[2] > 1 ? 1 : 0 ) ; 
    q(3,0) = 0.;     
    q(4,0) = 0.;    
    q(5,0) = 0.;


    bool EnoughShells = false;
    nshells = 0;
    Matrix<double> w;
    while(!EnoughShells) {
        nshells++;
        auto A_ = GradientMatrix(nshells, kmesh, ikshell); //computes A as in CPC 178 (2008) 685-599 between eq.25 and 26
        

        //compute w = A^{-1}*q
        auto invA = A_.pseudoinv();
        w = invA*q;

        //try A*w = qguess
        auto qguess = A_*w;

        //we have enough shells only when qguess is equal to q.
        EnoughShells = ( ( qguess - q ).norm() < threshold );
    }

    assert(w.get_nrows() == nshells);
    assert(w.get_ncols() == 1);

    Weight = mdarray<double, 1>({size_t(nshells)});
    for(int ishell = 0; ishell < nshells; ++ishell) {
        Weight(ishell) = w(ishell, 0);
    }
}


/*
 *         MAP
 *   j      0        1        2       3       4      5
 * alpha    x        y        z       x       x      y           
 * beta     x        y        z       y       z      z      
 */ 
int alpha( const size_t& j ) 
{
    if( j == 0 )       return 0;
    if( j == 1 )       return 1;
    if( j == 2 )       return 2;
    if( j == 3 )       return 0;
    if( j == 4 )       return 0;
    if( j == 5 )       return 1;
    return -1;
}

int beta( const size_t& j ) 
{
    if( j == 0 )       return 0;
    if( j == 1 )       return 1;
    if( j == 2 )       return 2;
    if( j == 3 )       return 1;
    if( j == 4 )       return 2;
    if( j == 5 )       return 2;
    return -1;
}

/*
 * This routine evalues A(j,s) = sum_{i\in s} b^{i}_\alpha b^{i}_\beta
 *
*/
Matrix<double> GradientMatrix(const size_t& nshells, const MeshGrid& kmesh, const std::vector<std::vector<int>>& ik_sorted)
{
    Matrix<double> A(6, nshells);
    A.fill(0.);

    for( size_t j=0; j<6; ++j ) { // j is a one d index for alpha, beta
        auto alpha_ = alpha(j);
        auto beta_  = beta(j);

        for( size_t ishell = 0; ishell < nshells; ishell++) {
            for ( auto& ikshell : ik_sorted[ishell] ) {// go through indices in shell
                auto& k = kmesh[ikshell].get("Cartesian");
                A(j, ishell) += k[alpha_]*k[beta_];
            }       
        }
    }
    return A;
}

std::vector<std::vector<std::vector<int>>> Find_kpb(const MeshGrid& kmesh, const std::vector<std::vector<int>>& ikshell)
{
    std::vector<std::vector<std::vector<int>>> ikpb(kmesh.get_TotalSize());

    std::cout << "Calculating map ik, ib ---> ikpb...\n";
    //find k+b in kmesh for every k, every b
    #pragma omp parallel for schedule(static)
    for( int ik = 0; ik < kmesh.get_TotalSize(); ++ik) {
        std::cout << "ik: " << ik << "/" << kmesh.get_TotalSize() << " in thread " << omp_get_thread_num() << "/" << omp_get_num_threads() << std::endl;
        ikpb[ik].resize(ikshell.size());
        for( int ishell = 0; ishell < ikshell.size(); ++ishell ) {
            ikpb[ik][ishell].resize( ikshell[ishell].size() );
            for( int ib = 0; ib < ikshell[ishell].size(); ++ib ) {
                ikpb[ik][ishell][ib] = kmesh.find(kmesh[ik] + kmesh[ikshell[ishell][ib]]);
            }
        }
    }
    std::cout << "DONE!";
    return ikpb;
}

