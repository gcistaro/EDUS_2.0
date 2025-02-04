template <typename Scalar_T>
void SumWithProduct(Operator<std::complex<double>>& Output_, 
                    const Scalar_T& FirstScalar_, 
                    const Operator<std::complex<double>>& FirstAddend_, 
                    const Scalar_T& SecondScalar_, 
                    const Operator<std::complex<double>>& SecondAddend_)
{
#ifdef EDUS_TIMERS
    PROFILE("SumWithProduct");
#endif
    auto& SpaceOfPropagation = Operator<std::complex<double>>::SpaceOfPropagation;

    auto& Output       = Output_      .get_Operator( SpaceOfPropagation );
    auto& FirstAddend  = FirstAddend_ .get_Operator( SpaceOfPropagation );
    auto& SecondAddend = SecondAddend_.get_Operator( SpaceOfPropagation );

    SumWithProduct(Output, FirstScalar_, FirstAddend, SecondScalar_, SecondAddend);
}