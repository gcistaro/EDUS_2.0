#include "MPIindex.hpp"



//constructor for splitting k point using fftw functions
template<size_t dim>
MPIindex<dim>::MPIindex( const std::array<int, 3>& ValuesToSplit__ )
{
    initialize(ValuesToSplit__ );
}

//for kpools
template<size_t dim>
void MPIindex<dim>::initialize( const std::array<int, 3>& ValuesToSplit__)
{
    //Global arrays are in range [0,values)
#ifdef NEGF_MPI    
    assert(ValuesToSplit__[0] > 1);//we dont want MPI to work without splitting anything, and fftw splits over this index.
    assert(ValuesToSplit__[1] > 1 || ValuesToSplit__[2] > 1 ); //because we know this fftw routines break for 1D FFT
                                                           //we can avoid this by making a choice on the fftw_mpi_local
                                                           //function to pick up
    std::ptrdiff_t dimensions[] = {std::ptrdiff_t( ValuesToSplit__[0] ), 
                                   std::ptrdiff_t( ValuesToSplit__[1] ),
                                   std::ptrdiff_t( ValuesToSplit__[2] )};
    std::ptrdiff_t local_0_start;
#endif
    ValuesToSplit[0] = ValuesToSplit__[0];
    ValuesToSplit[1] = ValuesToSplit__[1];
    ValuesToSplit[2] = ValuesToSplit__[2];

    std::cout << "Initializing multindex..\n";
    multindex.initialize(ValuesToSplit);
    //here we suppose the splitting happens with fftw. for the future must be fixed
    std::ptrdiff_t local_n0;

    GlobalRange_1D.first = 0;
    GlobalRange_1D.second = ValuesToSplit[0]*ValuesToSplit[1]*ValuesToSplit[2];

#ifdef NEGF_MPI
    auto alloc_local = fftw_mpi_local_size_many(int(dim), dimensions, 1,//howmany 
                                   FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD,
                                   &local_n0, &local_0_start);

    //initialize class variables using what we have obtained
    std::cout << "local_0_start " << local_0_start << " local_n0 " << local_n0 <<  std::endl;
    std::cout << "alloc_local  " << alloc_local <<  std::endl;
    LocalRange_1D.first  = multindex.oneDindex(int(local_0_start),0,0);
    LocalRange_1D.second = multindex.oneDindex(int(local_0_start + local_n0-1),ValuesToSplit[1]-1,ValuesToSplit[2]-1);

    RecommendedAllocate_fftw = alloc_local; //this is in general different than LocalRange.second-LocalRange.first
                                            //because of internals of fftw routines.
#else //trivial case, no splitting at all.
    LocalRange_1D = GlobalRange_1D;
    RecommendedAllocate_fftw = ValuesToSplit[0]*ValuesToSplit[1]*ValuesToSplit[2];
    local_n0 = ValuesToSplit[0];
#endif
    nlocal = local_n0*ValuesToSplit[1]*ValuesToSplit[2];
}

template<size_t dim>
inline std::ptrdiff_t MPIindex<dim>::glob1D_to_loc1D(const std::ptrdiff_t& global)
{
    assert( global >= GlobalRange_1D.first && global <= GlobalRange_1D.second );

    return global - LocalRange_1D.first;
}

template<size_t dim>
inline std::ptrdiff_t MPIindex<dim>::loc1D_to_glob1D(const std::ptrdiff_t& local)
{
    assert( local >= 0 && local < nlocal );
    return local + LocalRange_1D.first; 
}

template<size_t dim>
template<typename... Args>
inline std::ptrdiff_t MPIindex<dim>::globnD_to_glob1D(const Args&... args)
{
    return multindex.oneDindex(args...);
}

template<size_t dim>
inline std::array<int, dim> MPIindex<dim>::glob1D_to_globnD(const std::ptrdiff_t& glob1D)
{
    return multindex.nDindex(glob1D);
}

template<size_t dim>
template<typename... Args>
inline std::ptrdiff_t MPIindex<dim>::globnD_to_loc1D(const Args&... args)
{
    return multindex.oneDindex(args...) - LocalRange_1D.first;
}

template<size_t dim>
inline std::array<int, dim> MPIindex<dim>::loc1D_to_globnD(const std::ptrdiff_t& loc1D)
{
    assert( loc1D >= 0 && loc1D < nlocal );
    return multindex.nDindex( loc1D + LocalRange_1D.first );
}


template<size_t dim>
int MPIindex<dim>::get_RecommendedAllocate_fftw() const
{
    return RecommendedAllocate_fftw;
}

