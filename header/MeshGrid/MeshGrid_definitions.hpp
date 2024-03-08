template<Space space>
size_t MeshGrid<space>::counter_id = 0;
template<Space space>
std::map<std::array<size_t,3>, mdarray<int,2> > MeshGrid<space>::ConvolutionIndex1;
template<Space space>
std::map<std::array<size_t,3>, mdarray<int,2> > MeshGrid<space>::ConvolutionIndex2;

template<Space space>
MeshGrid<space>::MeshGrid(const MeshGrid<space>& mg)
{
    *this = mg;
}

template<Space space>
MeshGrid<space>& MeshGrid<space>::operator=(const MeshGrid<space>& mg)
{
    mesh.resize(mg.mesh.size());
    std::copy(&(mg.mesh[0]),&(mg.mesh[0])+mg.mesh.size(),&(mesh[0]));
    std::copy(&(mg.Size[0]),&(mg.Size[0])+3,&(Size[0]));
    std::copy(&(mg.type),&(mg.type)+1,&type);
    std::copy(&(mg.maxRadius), &(mg.maxRadius)+1, &(maxRadius));
    shellinfo.resize(mg.shellinfo.size());
    std::copy(mg.shellinfo.begin(), mg.shellinfo.end(), shellinfo.begin());
    EmptyVector = mg.EmptyVector;
    TotalSize = mg.TotalSize;
    id = mg.id;
    return *this;
}


template<Space space>
MeshGrid<space>::MeshGrid(MeshGrid<space>&& mg)
{
    *this = mg;
}

template<Space space>
MeshGrid<space>& MeshGrid<space>::operator=(MeshGrid<space>&& mg)
{
    mesh = std::move(mg.mesh);
    Size = std::move(mg.Size);
    type = std::move(mg.type);
    maxRadius = std::move(mg.maxRadius);
    shellinfo = std::move(mg.shellinfo);
    EmptyVector = std::move(EmptyVector);
    TotalSize = std::move(mg.TotalSize);
    id = std::move(mg.id);
    return *this;
}

template<Space space>
MeshGrid<space>::MeshGrid(const std::vector<Coordinate<space>>& ReadMesh) : type(read_), mesh(ReadMesh)
{
    TotalSize = mesh.size();
    ++counter_id;
    id = counter_id;
}

template<Space space>
MeshGrid<space>::MeshGrid(const std::array<int,3>& Size_)
{
    this->initialize(Size_);
}

template<Space space>
MeshGrid<space>::MeshGrid(const mdarray<double,2>& bare_mg, const std::string& KeyForBasis)
{
    this->initialize(bare_mg, KeyForBasis);
}


template<Space space>
MeshGrid<space>::MeshGrid(const double& Radius_, const std::array<double, 3>& resolution)
{
    this->initialize(Radius_, resolution);
}

template<Space space>
MeshGrid<space>::MeshGrid(const std::vector<Coordinate<space>>& PathPoint, const double& resolution)
{
    this->initialize(PathPoint, resolution);
}

template<Space space>
void MeshGrid<space>::initialize(const std::vector<Coordinate<space>>& PathPoint, const double& resolution)
{
    int NumberOfLines = PathPoint.size()-1;
    type = path;

    mesh.push_back(PathPoint[0]);
    for(int iline=0; iline<NumberOfLines; ++iline){
        int NumberOfPointsInLine = (PathPoint[iline+1]-PathPoint[iline]).norm()/resolution;
        for(int it=1; it<NumberOfPointsInLine; ++it){
            double t = double(it)/(NumberOfPointsInLine-1);
            mesh.push_back( (1-t)*PathPoint[iline] + t*PathPoint[iline+1] );
        }
    }
    TotalSize = mesh.size();
}

template<Space space>
void MeshGrid<space>::initialize(const std::array<int,3>& Size_)
{
    //cubic grids are always [0, Size-1] so that are fft friendly
    id = ++counter_id;

    type=cube;
    Size=Size_;
    TotalSize = Size[0]*Size[1]*Size[2];
    mesh.resize(TotalSize);
    //std::cout << "TotalSize: " << TotalSize << std::endl;
    //std::cout << "Size: " << Size[0] << " " << Size[1] << " " << Size[2]<<std::endl;
    auto KeyForBasis = "LatticeVectors";
    //defining vectors of the mesh
    switch(space){
        case(k):
        {
            for(int ix0=0; ix0<Size[0]; ix0++){
                for(int ix1=0; ix1<Size[1]; ix1++){
                    for(int ix2=0; ix2<Size[2]; ix2++){
                        int index = ix2 + Size[2] * (ix1 + Size[1] * ix0);
                        mesh[index].initialize(double(ix0)/Size[0],
                                               double(ix1)/Size[1],
                                               double(ix2)/Size[2],
                                               KeyForBasis);
                    }
                }
            }
            break;
        }
        case(R):
        {
            for(int ix0=0; ix0<Size[0]; ix0++){
                for(int ix1=0; ix1<Size[1]; ix1++){
                    for(int ix2=0; ix2<Size[2]; ix2++){
                        int index = ix2 + Size[2] * (ix1 + Size[1] * ix0);
                        mesh[index].initialize(double(ix0),
                                               double(ix1),
                                               double(ix2),
                                               KeyForBasis);
                    }
                }
            }
            break;
        }
    }
}

template<Space space>
void MeshGrid<space>::initialize(const mdarray<double,2>& bare_mg, const std::string& KeyForBasis)
{
    id = ++counter_id;
    type = read_;

    assert( bare_mg.get_Size(1) == 3 );
    for(int i=0; i<bare_mg.get_Size(0); ++i){
        mesh.push_back(Coordinate<space>(bare_mg(i,0), bare_mg(i,1), bare_mg(i,2), KeyForBasis));
    }
    TotalSize = mesh.size();
}

template<Space space>
void MeshGrid<space>::initialize(const double& Radius_, const std::array<double, 3>& resolution)
{
    PROFILE("MeshGrid::initializeSphere");
    //assert(space == R);
    id = ++counter_id;
    type = sphere;

    //defining vectors of the mesh
    double i0_threshold, i1_threshold, i2_threshold;

    i0_threshold = 50*resolution[0];
    i1_threshold = 50*resolution[1];
    i2_threshold = 50*resolution[2];
    for(int i0=-i0_threshold; i0<=i0_threshold; i0++){
        for(int i1=-i1_threshold; i1<=i1_threshold; i1++){
            for(int i2=-i2_threshold; i2<=i2_threshold; i2++){
                auto v = Coordinate<space>(double(i0)/resolution[0], double(i1)/resolution[1], double(i2)/resolution[2], "LatticeVectors");
                if(v.norm() < Radius_+threshold){
                    mesh.push_back(v);
                }
            }//end i2 loop
        }//end i1 loop
    }//end i0 loop

    std::sort(mesh.begin(), mesh.end(), [](const auto& v1, const auto& v2){ return v1.norm() < v2.norm();});
    TotalSize = mesh.size();

    maxRadius = Radius_;
    //here we order everything by shells
    double maxRadius = -1.;
    for(int iR=0; iR<mesh.size(); iR++){
        if(abs(mesh[iR].norm()-maxRadius) > 1.e-07){
            ShellInfo si;
            si.Radius = mesh[iR].norm();
            si.NumberOfPoints = 1;
            si.StartingIndex = iR;
            shellinfo.push_back(si);
            maxRadius = mesh[iR].norm();
        }
        else{
            ++shellinfo[shellinfo.size()-1].NumberOfPoints;
        }
    }
}


template<Space space>
std::array<double,3> MeshGrid<space>::get_absmax() const
{
    std::array< std::vector<double>, 3> absCoord;
    for (auto& absCoord_ : absCoord){
        absCoord_.reserve(mesh.size());
    }
    for(auto& vector : mesh){
        auto& CrystalCoordinate = vector.get("LatticeVectors");
        absCoord[0].push_back(abs(CrystalCoordinate[0]));
        absCoord[1].push_back(abs(CrystalCoordinate[1]));
        absCoord[2].push_back(abs(CrystalCoordinate[2]));
    }
    std::array<double, 3> maxCoord;
    for(int ix=0; ix<3; ++ix){
        maxCoord[ix] = *(std::max_element(absCoord[ix].begin(), absCoord[ix].end()));
    }
    return maxCoord;
}


template<Space space, Space space_>
MeshGrid<space_> fftPair(const MeshGrid<space>& KnownMG)
{
    PROFILE("MeshGrid::fftPair");
    assert(space_ != space);

    MeshGrid<space_> fftMeshGrid;

    switch(KnownMG.type)
    {
        case(sphere): case (read_):
        {
            auto maxCoord = KnownMG.get_absmax();
            std::array<int,3> maxCoordInt;
            for(int ix=0; ix<3; ix++){
                maxCoordInt[ix] = 2*int(maxCoord[ix])+1;
		        maxCoordInt[ix] += (maxCoordInt[ix]==0);
            }
            fftMeshGrid.initialize(maxCoordInt);
            break;
        }
        case(cube):
        {
            fftMeshGrid.initialize(KnownMG.get_Size());
            break;
        }
    }
    return fftMeshGrid;
}

template<Space space>
MeshGrid<space> Opposite(const MeshGrid<space>& MG)
{
    std::vector<Coordinate<space>> MGminus(MG.get_TotalSize());
    for(int i=0; i<MG.get_TotalSize(); i++){
        MGminus[i] = -MG[i];//MG.reduce(-MG[i]);
    }
    return MeshGrid<space>(MGminus);
}

template<Space space>
const Coordinate<space>& MeshGrid<space>::operator[](const int& i) const
{
    if(i < 0){
        return EmptyVector;
    }
    return mesh[i];
}

template<Space space>
Coordinate<space>& MeshGrid<space>::operator[](const int& i)
{
    return const_cast<Coordinate<space>&>(static_cast<const MeshGrid<space>&>(*this)[i]);
}


template<Space space>
int MeshGrid<space>::find(const Coordinate<space>& v) const
{
    int index = -1;
    //std::cout << "TYPE:: "<< type << std::endl;
    switch(type)
    {
        case cube:
        {
            auto v_reduced = reduce(v);
            auto notcart = v_reduced.get("LatticeVectors");
            index = int(round(notcart[2] + Size[2] * (notcart[1] + Size[1] * notcart[0])));//warning! round is important, if not it will get the integer part, usually wrong.
            break;
	    }
        case sphere:
        {
            if(v.norm() > maxRadius){
                return -1;
            }
            auto shellinfo_iterator = std::find_if(shellinfo.begin(), shellinfo.end(), [&v](const auto& si){ return abs(si.Radius-v.norm()) < 1.e-07; });
            auto BeginShell = mesh.begin()+(*shellinfo_iterator).StartingIndex;
            auto EndShell = mesh.begin()+(*shellinfo_iterator).StartingIndex+(*shellinfo_iterator).NumberOfPoints;
            auto v_iterator = std::find_if(BeginShell, EndShell, [&v](const auto& v_){return (v_-v).norm() < 1.e-07;});
            assert(v_iterator != EndShell);
            index = v_iterator-mesh.begin();
            break;
        }
        case read_:
        {

            auto v_iterator = std::find_if(mesh.begin(), mesh.end(), [&v](const auto& v_){return (v_-v).norm() < 1.e-07;});
            if(v_iterator == mesh.end()){
                return -1;
            }
            index = v_iterator-mesh.begin();
            break;
        }
    }
    if(index >= TotalSize) index = -1;
    return index;
}

template<Space space>
Coordinate<space> MeshGrid<space>::reduce(const Coordinate<space>& v) const
{
    assert(type == cube);
    auto notcart = v.get("LatticeVectors");
    std::array<int,3> grid_limit = (space==k) ? std::array<int,3>{1, 1, 1}
                                              : Size;

    for(int ix=0; ix<3; ix++){
        while(notcart[ix]<0 && abs(notcart[ix])>threshold){
            notcart[ix] += grid_limit[ix];
        }
    }

    Coordinate<space> v_reduced;
    v_reduced.initialize(std::fmod(notcart[0], grid_limit[0]),
                         std::fmod(notcart[1], grid_limit[1]),
                         std::fmod(notcart[2], grid_limit[2]),
                         "LatticeVectors");
    return v_reduced;
}

template<Space space>
const std::vector<Coordinate<space>>& MeshGrid<space>::get_mesh() const
{
    return mesh;
}


template<Space space>
const std::array<int,3>& MeshGrid<space>::get_Size() const
{
    return Size;
}

template<Space space>
size_t MeshGrid<space>::get_TotalSize() const
{
    return TotalSize;
}

template<Space space>
size_t MeshGrid<space>::get_id() const
{
    return id;
}


template<Space space>
void MeshGrid<space>::Calculate_ConvolutionIndex1(const MeshGrid& m1, const MeshGrid& m2, const MeshGrid& m3)
{
    PROFILE("MeshGrid::CalculateConvolutionIndex");
    auto i1 = m1.get_id();
    auto i2 = m2.get_id();
    auto i3 = m3.get_id();

    auto& ci = ConvolutionIndex1[{i1,i2,i3}];

    ci = mdarray<int,2>({m1.get_TotalSize(), m3.get_TotalSize()});
    std::cout << "Calculating convolution indices for " << i1 << " " << i2 << " " << i3 << std::endl;

    #pragma omp parallel for schedule(dynamic)
    for(int iR1=0; iR1<m1.get_TotalSize(); iR1++){
        if(!(iR1%200)){
            std::cout << "Calculate_ci::: iR1: " << iR1 << " / " << m1.get_TotalSize() << std::endl;
        }
        for(int iR3=0; iR3<m3.get_TotalSize(); iR3++){
            //std::cout << std::endl;
            //    std::cout << "iR1 << << iR3: " <<iR1 << " " << " " << iR3 << std::endl;
            //    std::cout << "m1[iR1]        : ";
            //    std::cout << m1[iR1].get("LatticeVectors");
            //    std::cout << "m3[iR3]        : ";
            //    std::cout << m3[iR3].get("LatticeVectors");
            //    std::cout << "m1[iR1]-m3[iR3]: ";
            //    std::cout << std::setprecision(15) << (m1[iR1]-m3[iR3]).get("LatticeVectors");

                ci(iR1, iR3) = m2.find(m1[iR1]-m3[iR3]);
            //auto& ciiR = ci(iR1,iR3);
            //    std::cout << "m2[ci]         : ";
            //if(ciiR >=0){
            //    std::cout << m2[ciiR].get("LatticeVectors");
            //    if((m1[iR1]-m3[iR3]-m2[ciiR]).norm() > threshold){
            //        std::cout << "HAVE A CONTROL OF THIS!\n";
            //    }
            //}
            //else{
            //    std::cout << "NOT FOUND\n";
            //}

            //    std::cout << iR1 << " " << iR3 << " ci(iR1,iR3): " << ci(iR1,iR3) << std::endl;
        }
    }
}


template<Space space>
void MeshGrid<space>::Calculate_ConvolutionIndex2(const MeshGrid& m1, const MeshGrid& m2, const MeshGrid& m3)
{

    PROFILE("MeshGrid::CalculateConvolutionIndex");
    auto i1 = m1.get_id();
    //std::cout << "Calculate_ConvolutionIndex " << i1<< std::endl;
    auto i2 = m2.get_id();
    //std::cout << "Calculate_ConvolutionIndex " << i2<< std::endl;
    auto i3 = m3.get_id();
    //std::cout << "Calculate_ConvolutionIndex " << i3<< std::endl;

    auto& ci = ConvolutionIndex2[{i1,i2,i3}];

    ci = mdarray<int,2>({m1.get_TotalSize(), m2.get_TotalSize()});
    std::cout << "Calculating convolution indices 2 for " << i1 << " " << i2 << " " << i3 << std::endl;
    //std::cout << "m1.get_TotalSize( ) " << m1.get_TotalSize() << " m3.get_TotalSize( ) " << m3.get_TotalSize() << std::endl;
    //std::cout << "Done.\n"<< std::endl;
    //The following openmp statement has been tested in one case.
    #pragma omp parallel for schedule(dynamic)
    for(int iR1=0; iR1<m1.get_TotalSize(); iR1++){
        if(!(iR1%200)){
            std::cout << "Calculate_ci::: iR1: " << iR1 << " / " << m1.get_TotalSize() << std::endl;
        }
        for(int iR2=0; iR2<m2.get_TotalSize(); iR2++){
        //    std::cout << std::endl;
        //        std::cout << "iR1 << << iR2: " <<iR1 << " " << " " << iR2 << std::endl;
        //        std::cout << "m1[iR1]        : ";
        //        std::cout << m1[iR1].get("LatticeVectors");
        //        std::cout << "m2[iR2]        : ";
        //        std::cout << m2[iR2].get("LatticeVectors");
        //        std::cout << "m1[iR1]+m2[iR2]: ";
        //        std::cout << std::setprecision(15) << (m1[iR1]+m2[iR2]).get("LatticeVectors");

                ci(iR1, iR2) = m3.find(m1[iR1]-m2[iR2]);//m3.find(m1[iR1]+m2[iR2]);
        //    auto& ciiR = ci(iR1,iR2);
        //        std::cout << "m2[ci]         : ";
        //    if(ciiR >=0){
        //        std::cout << m3[ciiR].get("LatticeVectors");
        //        if((m1[iR1]+m2[iR2]-m3[ciiR]).norm() > threshold){
        //            std::cout << "HAVE A CONTROL OF THIS!\n";
        //        }
        //    }
        //    else{
        //        std::cout << "NOT FOUND\n";
        //    }

            //    std::cout << iR1 << " " << iR3 << " ci(iR1,iR3): " << ci(iR1,iR3) << std::endl;
        }
    }}

template<Space space>
std::pair<int,int> MeshGrid<space>::get_shellindices(int shellNumber) const
{
    std::pair<int,int> indices;
    if(shellNumber < 0){
        shellNumber += shellinfo.size();
    }
    assert(shellNumber >= 0 && shellNumber < shellinfo.size());

    indices.first = shellinfo[shellNumber].StartingIndex;
    indices.second = indices.first + shellinfo[shellNumber].NumberOfPoints;

    return indices;

}
