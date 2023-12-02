template<Space space>
size_t MeshGrid<space>::counter_id = 0;
template<Space space>
std::map<std::array<size_t,3>, mdarray<int,2> > MeshGrid<space>::ConvolutionIndex;

template<Space space>
MeshGrid<space>& MeshGrid<space>::operator=(const MeshGrid<space>& mg)
{
    std::copy(&(mg.mesh[0]),&(mg.mesh[0])+mg.mesh.size(),&(mesh[0]));
    std::copy(&(mg.Size[0]),&(mg.Size[0])+3,&(Size[0]));
    std::copy(&(mg.type),&(mg.type)+1,&type);
    return *this;
}

template<Space space>
MeshGrid<space>& MeshGrid<space>::operator=(MeshGrid<space>&& mg)
{ 
     mesh = std::move(mg.mesh);
     Size = std::move(mg.Size);
     type = std::move(mg.type);
     return *this;
}

template<Space space>
MeshGrid<space>::MeshGrid(const std::vector<Coordinate<space>>& ReadMesh) : type(read), mesh(ReadMesh)
{
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
MeshGrid<space>::MeshGrid(const double& Radius_)
{
    this->initialize(Radius_);
}


template<Space space>
void MeshGrid<space>::initialize(const std::array<int,3>& Size_)
{
    id = ++counter_id;

    type=cube;
    Size=Size_;
    TotalSize = Size[0]*Size[1]*Size[2];
    mesh.resize(TotalSize);

    //defining vectors of the mesh
    auto BasisKey = "LatticeVectors";
    for(int ix0=0; ix0<Size[0]; ix0++){
        for(int ix1=0; ix1<Size[1]; ix1++){                   
            for(int ix2=0; ix2<Size[2]; ix2++){
                int index = ix2 + Size[2] * (ix1 + Size[1] * ix0);
                mesh[index].initialize(double(ix0),
                                       double(ix1),
                                       double(ix2),
                                       BasisKey);
            }
        }
    }
}

template<Space space>
void MeshGrid<space>::initialize(const mdarray<double,2>& bare_mg, const std::string& KeyForBasis)
{
    id = ++counter_id;
    type = read;

    assert( bare_mg.get_Size(1) == 3 );
    for(size_t i=0; i<bare_mg.get_Size(0); ++i){
        mesh.push_back(Coordinate<space>(bare_mg(i,0), bare_mg(i,1), bare_mg(i,2), KeyForBasis));
    }
    TotalSize = mesh.size();
}

/*
    in the following we use the inequality:
    ||i0 a0 + i1 a1 + i2 a2|| <= |i0| ||a0|| + |i1| ||a1|| + |i2| ||a2|| <=R
    so, starting from i0 we can get an upper boundary (overestimated in general) that asserts to get all the points in a sphere
    |i0max| = R/||a0||
    and for i1 and i2:
    |i1max| = (R-|i0| ||a0||)/||a1||
    |i2max| = (R-|i0| ||a0|| - |i1| ||a1||)/||a2||

*/
template<Space space>
void MeshGrid<space>::initialize(const double& Radius_)
{
    id = ++counter_id;
    type = sphere;
    

    //defining vectors of the mesh
    auto& MetricTensor_ = Coordinate<space>::get_Basis("LatticeVectors").get_MetricTensor();
    double i0_threshold, i1_threshold, i2_threshold;
    double fractPart;
    auto norm0 = std::sqrt(MetricTensor_(0,0));
    auto norm1 = std::sqrt(MetricTensor_(1,1));
    auto norm2 = std::sqrt(MetricTensor_(2,2));
    std::cout << "norm0 << norm1 << norm2 << " ;
    std::cout << norm0 << " " << norm1 << " " << norm2 << std::endl;

    fractPart = std::modf(Radius_/norm0, &i0_threshold);
    i0_threshold = abs(i0_threshold)+2; //this ensures to have a good upper bound

    fractPart = std::modf(Radius_/norm1, &i1_threshold);
    i1_threshold = abs(i1_threshold)+2; //this ensures to have a good upper bound

    fractPart = std::modf((Radius_+1.e-07)/norm2, &i2_threshold);
    i2_threshold = abs(i2_threshold)+2; //this ensures to have a good upper bound

    i0_threshold = 100;
    i1_threshold = 100;
    i2_threshold = 100;
    for(int i0=-i0_threshold; i0<=i0_threshold; i0++){
        for(int i1=-i1_threshold; i1<=i1_threshold; i1++){            
            for(int i2=-i2_threshold; i2<=i2_threshold; i2++){
                auto v = Coordinate<space>(double(i0), double(i1), double(i2), "LatticeVectors");
                if(v.norm() < Radius_+1.e-07){
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
    switch(type)
    {
        case cube:
        {
            auto v_reduced = reduce(v);
            auto notcart = v_reduced.get("LatticeVectors");
            index = int(notcart[2] + Size[2] * (notcart[1] + Size[1] * notcart[0]));
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
        case read:
        {
            auto v_iterator = std::find_if(mesh.begin(), mesh.end(), [&v](const auto& v_){return (v_-v).norm() < 1.e-07;});
            if(v_iterator == mesh.end()){
                return -1;
            }
            index = v_iterator-mesh.begin();
            break;
        }
    }
    return index;
}
        
template<Space space>
Coordinate<space> MeshGrid<space>::reduce(const Coordinate<space>& v) const
{
        auto notcart = v.get("LatticeVectors");
        for(int ix=0; ix<3; ix++){
            while(notcart[ix]<0 && abs(notcart[ix])>threshold){
                notcart[ix] += Size[ix];
            }
        }
        Coordinate<R> v_reduced;
        v_reduced.initialize(std::fmod(notcart[0],Size[0]),
                             std::fmod(notcart[1],Size[1]),
                             std::fmod(notcart[2],Size[2]),
                             "LatticeVectors");
        //std::cout << "I am inside reduce.. this is the vectors v_reduced at the end:\n";
        //std::cout << v_reduced;
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
int MeshGrid<space>::get_id() const
{
    return id;
}


template<Space space>
void MeshGrid<space>::Calculate_ConvolutionIndex(const MeshGrid& m1, const MeshGrid& m2, const MeshGrid& m3)
{
    auto i1 = m1.get_id();
    auto i2 = m2.get_id();
    auto i3 = m3.get_id();

    auto& ci = ConvolutionIndex[{i1,i2,i3}];
    ci = mdarray<int,2>({m1.get_TotalSize(), m3.get_TotalSize()});

    //The following openmp statement has been tested in one case.
    #pragma omp parallel for schedule(dynamic)
    for(int iR1=0; iR1<m1.get_TotalSize(); iR1++){
        std::cout << "Get convolution index ... " << iR1 << std::endl; 
        for(int iR3=0; iR3<m1.get_TotalSize(); iR3++){
                //std::cout << "iR1 << << iR3: " <<iR1 << " " << " " << iR3 << std::endl;
                //std::cout << m1[iR1].get("LatticeVectors") << m1[iR3].get("LatticeVectors");
                //std::cout << "m1[iR1]-m3[iR3]: ";
                //std::cout << std::setprecision(15) << (m1[iR1]-m1[iR3]).get("LatticeVectors") << std::endl;
                //std::cout << "m1[iR1]-m3[iR3].norm(): ";
                //std::cout << std::setprecision(15) << (m1[iR1]-m1[iR3]).norm() << std::endl;

                ci(iR1, iR3) = m2.find(m1[iR1]-m3[iR3]);
                
                //std::cout << "ci(iR1,iR3): " << ci(iR1,iR3) << std::endl;
        }
    }
}
