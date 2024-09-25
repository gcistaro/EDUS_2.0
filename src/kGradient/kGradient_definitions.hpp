
/*
   This function calculates the gradient supposing you have a DerivativeFunction already allocated!
   in k space, we use the formula from Mostofi2008; in R space, we calculate it as:
   Derivative(R) = i*R*Function(R)
*/
template<typename T, typename U>
void kGradient::Calculate(T& DerivativeFunction, const T& Function, 
                          const U& direction, const bool& EraseOutput) const
{
    PROFILE("kGradient::Calculate");

    if( EraseOutput ) {
        DerivativeFunction.fill(0.);
    }

    /* case 1-> we calculate the gradient directly in k */
    if( kmesh && !Rmesh ) {
        #pragma omp parallel for schedule(static)
        for( int ik = 0; ik < mpindex.nlocal; ++ik ) {
            for( int ishell = 0; ishell < Weight.get_Size(0); ++ishell ) {
                for( int ib = 0; ib < int( ikshell[ishell].size() ); ++ib ) {
                    auto& Bvector = (*kmesh)[ikshell[ishell][ib]];
                    auto bdotu = Bvector.dot(direction);
                    auto& ikpb_ = ikpb[ik][ishell][ib];
                    DerivativeFunction[ik] += Weight(ishell)*bdotu*Function[ikpb_];  
                }
            }
        }
    }
    /* case 2-> we calculate the gradient in R. It is up to the user to make sure Function
       is defined in R and to go back to k after doing it. TODO!! Change this!!
    */
    else if ( Rmesh && !kmesh ) {
        #pragma omp parallel for schedule(static)
        for( int iR_loc = 0; iR_loc < mpindex.nlocal; ++iR_loc ) {
            auto iR_global = (static_cast<MPIindex<3>>(mpindex)).loc1D_to_glob1D(iR_loc);
            DerivativeFunction[iR_loc] += im*(*Rmesh)[iR_global].dot(direction)*Function[iR_loc];
        }
    }
}