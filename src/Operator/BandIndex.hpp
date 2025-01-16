
#ifndef BANDINDEX_HPP
#define BANDINDEX_HPP

/*********************************************************************************************
 *   numbering:
 *   [0   1   2    3 ]
 *   [.   4   5    6 ]
 *   [.   .   7    8 ]
 *   [.   .   .    9 ]
 *
 *   ......
 *   Total number of elements before row m:
 *   b + (b-1) + ... + (b-m-1) = \sum_{i=0}^{m-1} (b-i) = m*b-(m-1)*(m)/2  .... VALID FOR m>1
 *   This is the starting index of row m.
 *
 *   The final index is
 *   StartingIndex + #elements in row = m*b-(m-1)*m/2 + (b-m-1) 
 *
*/ 

class BandIndex
{
    private:
        size_t NumberOfBands;
        size_t oneDNumberOfBands;

        std::vector< std::pair<size_t,size_t> > RowIndexBoundary;

    public:
    size_t StartingIndex(const auto& row_)
    {
        size_t StartingIndex = (row_ == 0) ? 0 
                                           :row_*NumberOfBands - row_*(row_-1)/2;
        return StartingIndex;
    };

    size_t RowIndex(const size_t& oneDband)
    {
        auto RowIndexIterator_ = std::find_if(RowIndexBoundary.begin(), RowIndexBoundary.end(),
                                             [&](const auto& RowPair){return oneDband >= RowPair.first && 
                                                                             oneDband <= RowPair.second;
                                                                      });
        return RowIndexIterator_ - RowIndexBoundary.begin();
    };

    BandIndex(){};
    BandIndex(const size_t& NumberOfBands__)
    {
        initialize(NumberOfBands__);
    };

    void initialize(const size_t& NumberOfBands__)
    {
        NumberOfBands = NumberOfBands__;
        oneDNumberOfBands = StartingIndex(NumberOfBands);
        RowIndexBoundary.resize(NumberOfBands);
        for(int i=0; i<NumberOfBands; ++i){
            RowIndexBoundary[i].first = StartingIndex(i);
            RowIndexBoundary[i].second = StartingIndex(i) + (NumberOfBands - i) - 1;
        }

    }
    size_t oneDband(const size_t& bnd1, const size_t& bnd2)
    {
        assert(bnd1 <= bnd2);
        assert(bnd1 <= NumberOfBands && bnd2 <= NumberOfBands);

        return RowIndexBoundary[bnd1].first + (bnd2-bnd1);
    };

    std::pair<size_t, size_t> twoDband(const size_t& oneDband__)
    {
        assert(oneDband__ < oneDNumberOfBands);
        std::pair<size_t, size_t> twoDband_;
        twoDband_.first = RowIndex(oneDband__);
        twoDband_.second = oneDband__ - RowIndexBoundary[twoDband_.first].first + twoDband_.first;
        return twoDband_;
    };

    size_t& get_oneDNumberOfBands(){ return oneDNumberOfBands; };
};

#endif