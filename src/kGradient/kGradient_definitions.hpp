
template<typename T, typename U>
void kGradient::Calculate(T& DerivativeFunction, const T& Function, 
                          const U& direction, const bool& EraseOutput) const
{
    if( EraseOutput ) {
        DerivativeFunction.fill(0.);
    }


    #pragma omp parallel for schedule(static)
    for( int ik = 0; ik < kmesh->get_TotalSize(); ++ik ) {
        for( int ishell = 0; ishell < Weight.get_Size(0); ++ishell ) {
            for( int ib = 0; ib < ikshell[ishell].size(); ++ib ) {
                auto& Bvector = (*kmesh)[ikshell[ishell][ib]];
                auto bdotu = Bvector.dot(direction);
                if( std::abs( bdotu ) > threshold ) {
                    auto& ikpb_ = ikpb[ik][ishell][ib];
                    DerivativeFunction[ik] += Weight(ishell)*bdotu*Function[ikpb_];  
                }
            }
        }
    }
}