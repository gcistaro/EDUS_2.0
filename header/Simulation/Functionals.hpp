std::function<void(BlockMatrix<std::complex<double>, R>&)> InitialCondition = [&](BlockMatrix<std::complex<double>, R>& DM)
{
    PROFILE("RK::InitialCondition");
    //DM.Operator_k.initialize(Uk.get_nblocks(), Uk.get_nrows(), Uk.get_ncols());
    DensityMatrix.Operator_k.fill(0.);
    //filling matrix in Bloch gauge
    for(int ik=0; ik<Uk.get_nblocks(); ++ik){
        for(int iband=0; iband<Uk.get_nrows(); iband++){
            if(this->Band_energies[ik](iband) < FermiEnergy-threshold){
                DensityMatrix.Operator_k(ik, iband, iband) = 1.;
            }
        }
    }
    DensityMatrix.lock_gauge(bloch);
    DensityMatrix.lock_space(k);
    DensityMatrix.go_to_wannier();
    DensityMatrix.go_to_R();


};

std::function<void(BlockMatrix<std::complex<double>, k>&)> InitialCondition_k = [&](BlockMatrix<std::complex<double>, k>& DM)
{
    PROFILE("RK::InitialCondition");
    //DM.Operator_k.initialize(Uk.get_nblocks(), Uk.get_nrows(), Uk.get_ncols());
    DensityMatrix.Operator_k.fill(0.);
    //filling matrix in Bloch gauge
    for(int ik=0; ik<Uk.get_nblocks(); ++ik){
        for(int iband=0; iband<Uk.get_nrows(); iband++){
            if(this->Band_energies[ik](iband) < FermiEnergy-threshold){
                DensityMatrix.Operator_k(ik, iband, iband) = 1.;
            }
        }
    }
    DensityMatrix.lock_gauge(bloch);
    DensityMatrix.lock_space(k);
    DensityMatrix.go_to_wannier();
    //DensityMatrix.go_to_R();
};

std::function<void(BlockMatrix<std::complex<double>, R>&, double const&, BlockMatrix<std::complex<double>, R> const&)> SourceTerm = 
[&](BlockMatrix<std::complex<double>, R>& Output, const double& time, const BlockMatrix<std::complex<double>, R>& Input)
{
    PROFILE("RK::SourceTerm");
    //Von Neumann Equations: idp/dt = [H, p]

    //DensityMatrix.go_to_k();
    //DensityMatrix.go_to_R();
    //DensityMatrix.go_to_k();
    //BlockMatrix<std::complex<double>, k> Out = DensityMatrix.get_Operator_k();
    //BlockMatrix<std::complex<double>, k> Out1 = DensityMatrix.get_Operator_k();
    //BlockMatrix<std::complex<double>, k> Out2 = DensityMatrix.get_Operator_k();
    //Out1.fill(0.);
    //Out2.fill(0.);
    //
//
    //multiply(Out1, -im, material.H.get_Operator_k(), DensityMatrix.get_Operator_k());
    //multiply(Out2, +im, DensityMatrix.get_Operator_k(), material.H.get_Operator_k());
    //for(int ik=0; ik<Out1.get_MeshGrid()->get_TotalSize(); ik++){
    //    SumWithProduct(Out[ik], 1., Out1[ik], 1., Out2[ik]); 
    //}
    //std::cout << std::endl << "Maximum of out: "<< " " << max(Out).data() << " " <<  max(Out)-Out.begin() << " " << *max(Out) << Out.Values[max(Out)-Out.begin()] << std::endl;    



    //Calculate_TDHamiltonian(time);
    //convolution1(Output, -im, H.get_Operator_R(), Input);
    //convolution1(Output, +im, Input, H.get_Operator_R() );
    Output.fill(0.*im);
    //std::cout <<"\n\n\n\n\n\n\n\n" << Input << "\n\n\n\n\n\n\n\n" ;
    convolution1(Output, -im, material.H.get_Operator_R(), Input);
    convolution1(Output, +im, Input, material.H.get_Operator_R() );
    
    std::cout << "Max of Input: " << *max(Input) << std::endl;
    auto index = int(max(Output)-Output.begin());
    std::cout << "Max of Output: " << *(max(Output)) << index << std::endl;

    //for(int iR=0; iR<Output.get_MeshGrid()->get_TotalSize(); ++iR){
    //    if(std::abs(max(Output[iR]))>1.e-05){
    //        std::cout << (*Output.get_MeshGrid())[iR].get("LatticeVectors")  << Output[iR] << std::endl;
    //    }
    //}

    //exit(0);

    //auto& m1 = *(Output.get_MeshGrid());
    //auto& m2 = *(const_cast<BlockMatrix<std::complex<double>,R>&>(Input).get_MeshGrid());
    //auto& m3 = *(material.H.get_Operator_R().get_MeshGrid());
    //auto& c1 = MeshGrid<R>::ConvolutionIndex1[{m1.get_id(), m3.get_id(), m2.get_id()}];
    //auto& c2 = MeshGrid<R>::ConvolutionIndex2[{m1.get_id(), m2.get_id(), m3.get_id()}];



    //for(int iR1=0; iR1<m1.get_mesh().size(); iR1++){
    //    if(std::abs(max(Output[iR1]))>threshold){
    //        std::cout << iR1 << " " << m1[iR1].get("LatticeVectors") << Output[iR1] << std::endl;
    //        for(int iR2=0; iR2<m2.get_mesh().size(); iR2++){
    //            if(c1(iR1,iR2) != -1){
    //                std::cout << "m2[iR2] " << m2[iR2].get("LatticeVectors") ;
    //                std::cout << Input[iR2] << std::endl;
    //            }
    //        }
    //        for(int iR2=0; iR2<m2.get_mesh().size(); iR2++){
    //            if(c1(iR1,iR2) != -1){
    //                std::cout << "TYPE1\n";
    //                std::cout << "m1[iR1] " << m1[iR1].get("LatticeVectors") ;
    //                std::cout << "m2[iR2] " << m2[iR2].get("LatticeVectors") ;
    //                std::cout << "m3[ci] " << m3[c1(iR1,iR2)].get("LatticeVectors") ;
    //                std::cout << material.H.get_Operator_R()[c1(iR1,iR2)]*Input[iR2] << std::endl;
    //                
    //            //std::cout << "m2[iR2] " << m2[iR2].get("LatticeVectors");
    //            //std::cout << "m3[c1] " << m3[c1(index, iR2)].get("LatticeVectors") << std::endl;
    //            }
    //            if(c2(iR1,iR2) != -1){
    //                std::cout << "TYPE2\n";
    //                std::cout << "m1[iR1] " << m1[iR1].get("LatticeVectors") ;
    //                std::cout << "m2[iR2] " << m2[iR2].get("LatticeVectors") ;
    //                std::cout << "m3[ci] " << m3[c2(iR1,iR2)].get("LatticeVectors") ;
    //                std::cout << Input[iR2]*material.H.get_Operator_R()[c2(iR1,iR2)];
    //            }
    //        }
    //            std::cout << std::endl << std::endl;
    //            std::cout << std::endl << std::endl;
    //    }
    //exit(0);
    //}
//    for(int iR2=0; iR2<m2.get_mesh().size(); iR2++){
//        if(c1(index,iR2) != -1){
//            std::cout << "m1[index] " << m1[index].get("LatticeVectors") ;
//            std::cout << "m2[iR2] " << m2[iR2].get("LatticeVectors");
//            std::cout << "m3[c1] " << m3[c1(index, iR2)].get("LatticeVectors") << std::endl;
//        }
//    }
    //std::cout << "Input: " << Input[0];
    //std::cout << "Output: " << Output[0];
    //auto& Rmesh = Output.get_MeshGrid()->get_mesh();
    //for(auto& R : Rmesh){
    //    std::cout << R.get("LatticeVectors")[0] << " " <<R.get("LatticeVectors")[1] << " " << R.get("LatticeVectors")[2]<< std::endl;
    //}
};


std::function<void(BlockMatrix<std::complex<double>, k>&, double const&, BlockMatrix<std::complex<double>, k> const&)> SourceTerm_k = 
[&](BlockMatrix<std::complex<double>, k>& Output, const double& time, const BlockMatrix<std::complex<double>, k>& Input)
{
    PROFILE("RK::SourceTerm");
    //Von Neumann Equations: idp/dt = [H, p]
    Output.fill(0.*im);
    Calculate_TDHamiltonian(time);
    multiply(Output, -im, material.H.get_Operator_k(), Input);
    //material.H.go_to_bloch();
    //std::cout << "materialH: " << material.H.get_Operator_k();
    
    multiply(Output, +im, Input, material.H.get_Operator_k(), 1.+0.*im );
    
    std::cout << "Max of Input: " << *(max(Input)) << std::endl;
    std::cout << "Max of Output: " << *(max(Output));// << int(max(Output)-Output.begin()) << std::endl;
    std::cout << "materialH: " << material.H.get_Operator_k();
    std::cout << "Input: " << Input;
    std::cout << "Output: " << Output;
};
