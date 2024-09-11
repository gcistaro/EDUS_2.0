#include "Simulation/Simulation.hpp"

bool Simulation::PrintObservables(const double& time) const
{
    return ( int( round( time/RK_object.get_ResolutionTime() ) ) % PrintResolution == 0 );
}

void Simulation::SettingUp_EigenSystem()
{
    PROFILE("Simulation::SetEigensystem");

    //------------------------go to H(k)--------------------------------------------
    auto& MasterkGrid = DensityMatrix.get_Operator_k().get_MeshGrid()->get_mesh();
    material.H.dft(MasterkGrid, +1);
    //------------------------------------------------------------------------------
    std::stringstream ss;
    ss << "rank_"<< mpi::Communicator::world().rank();
    std::ofstream processor_H(ss.str());
    for(int ik_loc=0; ik_loc< material.H.get_Operator_k().get_nblocks(); ++ik_loc) {
        processor_H << material.H.get_Operator_k()[ik_loc] << std::endl;
    }
    //--------------------solve eigen problem---------------------------------------
    auto& Uk = Operator<std::complex<double>>::EigenVectors;
    auto& UkDagger = Operator<std::complex<double>>::EigenVectors_dagger;

    material.H.get_Operator_k().diagonalize(Band_energies, Uk);
    //------------------------------------------------------------------------------

    //-------------------get U dagger-----------------------------------------------
    UkDagger.initialize(k, Uk.get_nblocks(), Uk.get_ncols(), Uk.get_nrows());
    for(int ik=0; ik<UkDagger.get_nblocks(); ++ik){
        for(int ir=0; ir<UkDagger.get_nrows(); ++ir){
            for(int ic=0; ic<UkDagger.get_ncols(); ++ic){
                UkDagger[ik](ir, ic) = conj(Uk[ik](ic,ir));
            }
        }
    }
    //------------------------------------------------------------------------------
};   


void Simulation::Calculate_TDHamiltonian(const double& time)
{
    PROFILE("SourceTerm::Calculate_TDHamiltonian");
    //--------------------get aliases for nested variables--------------------------------
    auto& H_ = H.get_Operator(SpaceOfPropagation);
    auto& H0_ = material.H.get_Operator(SpaceOfPropagation);
    auto& x_ = material.r[0].get_Operator(SpaceOfPropagation);
    auto& y_ = material.r[1].get_Operator(SpaceOfPropagation);
    auto& z_ = material.r[2].get_Operator(SpaceOfPropagation);
    auto las = laser(time).get("Cartesian");
    
    auto& ci = MeshGrid::ConvolutionIndex[{H0_.get_MeshGrid()->get_id(), 
                                              H_.get_MeshGrid()->get_id(), 
                                              Operator<std::complex<double>>::MeshGrid_Null->get_id()}];
    //-----------------------------------------------------------------------------------

    //--------------------------do initializations---------------------------------------
    H_.fill(0.);

    if(ci.get_Size(0) == 0 && SpaceOfPropagation == R) {
        MeshGrid::Calculate_ConvolutionIndex( *(H0_.get_MeshGrid()), 
                                                  *(H_.get_MeshGrid()), 
                                                  *(Operator<std::complex<double>>::MeshGrid_Null));
    }
    //---------------------------------------------------------------------------------

    //------------------------H(R) = H0(R) + E.r(R)-------------------------------------
    #pragma omp parallel for schedule(static) collapse(3)
    for(int iblock=0; iblock<H0_.get_nblocks(); ++iblock){
        for(int irow=0; irow<H0_.get_nrows(); ++irow){
            for(int icol=0; icol<H0_.get_ncols(); ++icol){
                //auto Hblock = ( SpaceOfPropagation == k ? iblock : ci(iblock, 0) );
                //H_(Hblock, irow, icol) = H0_(iblock, irow, icol)
                H_(iblock, irow, icol) = H0_(icol+2*irow+2*iblock)//iblock, irow, icol)
                                       + las[0]*x_(iblock, irow, icol)
                                       + las[1]*y_(iblock, irow, icol)
                                       + las[2]*z_(iblock, irow, icol);
                //std::cout << iblock << " " << irow << " " << icol << " " << H_(Hblock, irow, icol)-H0_(iblock, irow, icol) << std::endl;
            }
        }
    }
    //---------------------------------------------------------------------------------
}

void Simulation::Propagate()
{
    PROFILE("Simulation::Propagate");
    auto CurrentTime = RK_object.get_CurrentTime();
    //------------------------Print population-------------------------------------
    if(SpaceOfPropagation == R) DensityMatrix.go_to_k();
    //DensityMatrix.go_to_bloch();
    
    /*
    auto Population = TraceK(DensityMatrix.get_Operator_k());
    for(int ibnd=0; ibnd<Population.size(); ibnd++){
        os_Pop << std::setw(20) << std::setprecision(10) << Population[ibnd].real();
        os_Pop << ' ';
    }
    os_Pop << std::endl;

*/
    if( PrintObservables( CurrentTime) ){
        //print time 
        os_Time << CurrentTime << std::endl;
         //print population
/*
        std::vector<std::complex<double>> Population(DensityMatrix.get_Operator_k().get_nrows(), 0.);
        for(int ik=0; ik<DensityMatrix.get_Operator_k().get_nblocks(); ik++) {
            for(int ibnd=0; ibnd<Population.size(); ibnd++) {
                Population[ibnd] += DensityMatrix.get_Operator_k()[ik](ibnd,ibnd);
            }
        }
        Population[0] = 1.-Population[0]/double(DensityMatrix.get_Operator_k().get_MeshGrid()->get_TotalSize());
        Population[1] /= double(DensityMatrix.get_Operator_k().get_MeshGrid()->get_TotalSize());


        for(int ibnd=0; ibnd<Population.size(); ibnd++){
            os_Pop << std::setw(20) << std::setprecision(10) << Population[ibnd].real();
            os_Pop << ' ';
        }
        os_Pop << std::endl;
*/
        mdarray<std::complex<double>,2> Population;
        DensityMatrix.go_to_bloch();
        Population.initialize({DensityMatrix.get_Operator_k().get_nrows(), DensityMatrix.get_Operator_k().get_ncols()});
        Population.fill(0.);
        for(int ik=0; ik<DensityMatrix.get_Operator_k().get_nblocks(); ik++) {
            for(int irow=0; irow <DensityMatrix.get_Operator_k().get_nrows(); ++irow ) {
                for(int icol=0; icol <DensityMatrix.get_Operator_k().get_ncols(); ++icol ) {
                    Population(irow, icol) += DensityMatrix.get_Operator_k()[ik](irow, icol);
                }
            }
        }
        for(int irow=0; irow <DensityMatrix.get_Operator_k().get_nrows(); ++irow ) {
            for(int icol=0; icol <DensityMatrix.get_Operator_k().get_ncols(); ++icol ) {
                if(irow == 0 && icol == 0) Population(irow, icol) = double(DensityMatrix.get_Operator_k().get_nblocks())-Population(irow, icol);
                os_Pop << std::scientific <<  Population(irow, icol)/double(DensityMatrix.get_Operator_k().get_nblocks());
                os_Pop << ' ';
            }
        }
        os_Pop << std::endl;
        //print laser
        os_Laser << laser(RK_object.get_CurrentTime()).get("Cartesian");
	    //print stuff

        DensityMatrix.go_to_wannier();
    }

    //------------------------------------------------------------------------------


    //
    //std::array<std::complex<double>, 3> v;
    //for(auto ix : {0, 1, 2}){
    //    v[ix] = Trace(Velocity[ix].get_Operator_R(), DensityMatrix.get_Operator_R());
    //}
    //os_Pos << std::setw(20) << std::setprecision(8) << v[0].real();
    //os_Pos << std::setw(20) << std::setprecision(8) << v[0].imag();
    //os_Pos << std::setw(20) << std::setprecision(8) << v[1].real();
    //os_Pos << std::setw(20) << std::setprecision(8) << v[1].imag();
    //os_Pos << std::setw(20) << std::setprecision(8) << v[2].real();
    //os_Pos << std::setw(20) << std::setprecision(8) << v[2].imag();
    //os_Pos << std::endl;

    //DensityMatrix.go_to_wannier();
    //DensityMatrix.go_to_R();
    RK_object.Propagate();
}

void Simulation::Calculate_Velocity()
{
    for(int ix : {0, 1, 2}){
        Velocity[ix].get_Operator_R() = DensityMatrix.get_Operator_R();
    
        Velocity[ix].get_Operator_R().fill(0.);
        commutator(Velocity[ix].get_Operator_R(), -im, material.r[ix].get_Operator_R(), material.H.get_Operator_R());
        //part with R
        auto& ci = MeshGrid::ConvolutionIndex[{Velocity[ix].get_Operator_R().get_MeshGrid()->get_id(), 
                                              material.H.get_Operator_R().get_MeshGrid()->get_id(), 
                                              Operator<std::complex<double>>::MeshGrid_Null->get_id()}];
        if(ci.get_Size(0) == 0 ){
            MeshGrid::Calculate_ConvolutionIndex(*(Velocity[ix].get_Operator_R().get_MeshGrid()), 
                                                     *(material.H.get_Operator_R().get_MeshGrid()), 
                                                     *(Operator<std::complex<double>>::MeshGrid_Null));
        }

        for(int iblock=0; iblock<Velocity[ix].get_Operator_R().get_nblocks(); ++iblock){
            auto& R =(*(Velocity[ix].get_Operator_R().get_MeshGrid()))[iblock].get("Cartesian"); 
            for(int irow=0; irow<Velocity[ix].get_Operator_R().get_nrows(); ++irow){
                for(int icol=0; icol<Velocity[ix].get_Operator_R().get_ncols(); ++icol){
                    if(ci(iblock, 0) != -1){
                        Velocity[ix].get_Operator_R()(iblock, irow, icol) += -im*R[ix]*material.H.get_Operator_R()[ci(iblock,0)](irow, icol);
                    }
                }
            }
        }
    }
}


void Simulation::print_recap()
{
    std::cout << "*************  GRIDS: *************\n";
    std::cout << "Total number of grids:  " << MeshGrid::get_counter_id() << std::endl;
    std::cout << "              ";
    std::cout << "| id |";
    std::cout << " size  |\n";
    std::cout << "material grid ";
    std::cout << "|" << std::setw(4) << material.H.get_Operator_R().get_MeshGrid()->get_id() << "|";
    std::cout<<  std::setw(7)  <<  material.H.get_Operator_R().get_MeshGrid()->get_mesh().size() << "|"<< std::endl;
    std::cout << "DM grid (R)   ";
    std::cout << "|" << std::setw(4) << DensityMatrix.get_Operator_R().get_MeshGrid()->get_id() << "|";
    std::cout <<  std::setw(7) << DensityMatrix.get_Operator_R().get_MeshGrid()->get_mesh().size() << "|"<< std::endl;
    std::cout << "DM grid (k)   ";
    std::cout << "|" << std::setw(4) << DensityMatrix.get_Operator_k().get_MeshGrid()->get_id() << "|";
    std::cout <<  std::setw(7) << DensityMatrix.get_Operator_k().get_MeshGrid()->get_mesh().size() << "|"<< std::endl;
    laser.print_info();
}






void Simulation::print_grids()
{
    auto& MasterRgrid = DensityMatrix.get_Operator_R().get_MeshGrid();
    auto& kgrid = DensityMatrix.get_Operator_k().get_MeshGrid();
    auto Rfft = std::make_shared<MeshGrid>(fftPair(*kgrid));//Rspace as fourier space of kgrid
    std::ofstream os;
    os.open("MasterR.txt");
    for (auto &R : MasterRgrid->get_mesh())
    {
        os << R.get("Cartesian");
    }
    os.close();
    os.open("k.txt");
    for (auto &k : kgrid->get_mesh())
    {
        os << k.get("Cartesian");
    }
    os.close();
    auto MeshRfft = fftPair(*kgrid);
    os.open("R.txt");
    for (auto &R : Rfft->get_mesh())
    {
        os << R.get("Cartesian");
    }
    os.close();
}