
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
        for( int ik = 0; ik < kmesh->get_TotalSize(); ++ik ) {
            for( int ishell = 0; ishell < Weight.get_Size(0); ++ishell ) {
                for( int ib = 0; ib < ikshell[ishell].size(); ++ib ) {
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
        for( int iR = 0; iR < Rmesh->get_TotalSize(); ++iR ) {
            DerivativeFunction[iR] = im*2.*pi*(*Rmesh)[iR].dot(direction)*Function[iR];
        }
    }
}