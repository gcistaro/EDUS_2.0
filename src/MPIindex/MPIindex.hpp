
#ifndef MPIINDEX_HPP
#define MPIINDEX_HPP

#include "mpi.h"
#include "MultiIndex/MultiIndex.hpp"

template<size_t dim>
class MPIindex
{
    private:
        std::pair<std::ptrdiff_t, std::ptrdiff_t> GlobalRange_1D;
        std::pair<std::ptrdiff_t, std::ptrdiff_t> LocalRange_1D;
        int RecommendedAllocate_fftw;
        std::array<size_t, dim> ValuesToSplit;
        std::array<int, dim> Offset;
        MultiIndex<dim> multindex;

    public:
        MPIindex(){};

        MPIindex( const std::array<size_t, 3>& ValuesToSplit__);

        void initialize( const std::array<size_t, 3>& ValuesToSplit__ );
        
        inline std::ptrdiff_t glob1D_to_loc1D(const std::ptrdiff_t& global);
        inline std::ptrdiff_t loc1D_to_glob1D(const std::ptrdiff_t& local);

        template<typename... Args>
        inline std::ptrdiff_t globnD_to_glob1D(const Args&... args);

        template<typename... Args>
        inline std::ptrdiff_t globnD_to_loc1D(const Args&... args);

        inline std::array<int, dim> loc1D_to_globnD(const std::ptrdiff_t& loc1D);

        inline std::array<int, dim> glob1D_to_globnD(const std::ptrdiff_t& glob1D);

        inline int get_RecommendedAllocate_fftw() const;

        inline const auto& get_GlobalRange() const {return GlobalRange_1D; };
        inline auto& get_GlobalRange() {return GlobalRange_1D; };


        inline const auto& get_LocalRange() const {return LocalRange_1D; };
        inline auto& get_LocalRange() {return LocalRange_1D; };

};


#include "MPIindex/MPIindex_definitions.hpp"
#endif