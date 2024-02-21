template<typename T, size_t dim> 
mdarray<T,dim>::mdarray(const mdarray<T,dim>& ToBeCopied)
{
    *this = ToBeCopied;
}

//copy assignment 
template<typename T, size_t dim> 
mdarray<T,dim>& mdarray<T,dim>::operator=(const mdarray<T,dim>& ToBeCopied)
{
    //this->Ptr = std::make_unique<T[]>(ToBeCopied.TotalSize);
    (*this).~mdarray<T,dim>();
    this->Ptr = new T[ToBeCopied.TotalSize];
    this->Size = ToBeCopied.Size;
    this->TotalSize = ToBeCopied.TotalSize;
    this->Offset = ToBeCopied.Offset;
    std::copy(ToBeCopied.begin(), ToBeCopied.end(), this->begin());
    return *this;
}

//move constructor
template<typename T, size_t dim> 
mdarray<T,dim>::mdarray(mdarray<T,dim>&& ToBeMoved)
{
    (*this).~mdarray<T,dim>();
    this->Ptr = std::move(ToBeMoved.Ptr);
    this->Size = std::move(ToBeMoved.Size);
    this->TotalSize = std::move(ToBeMoved.TotalSize);
    if(this->TotalSize>0)    
    this->Offset = std::move(ToBeMoved.Offset);
    ToBeMoved.Ptr=nullptr;
    NotDestruct = ToBeMoved.NotDestruct;
}


//move assignment 
template<typename T, size_t dim> 
mdarray<T,dim>& mdarray<T,dim>::operator=(mdarray<T,dim>&& ToBeMoved)
{
    (*this).~mdarray<T,dim>();
    this->Ptr = std::move(ToBeMoved.Ptr);
    this->Size = std::move(ToBeMoved.Size);
    this->TotalSize = std::move(ToBeMoved.TotalSize);
    if(this->TotalSize>0)    
    this->Offset = std::move(ToBeMoved.Offset);
    ToBeMoved.Ptr=nullptr;
    NotDestruct = ToBeMoved.NotDestruct;
    return *this;
}


template<typename T, size_t dim> 
mdarray<T,dim>::mdarray(const std::array<size_t,dim>& Size_) 
{
    this->initialize(Size_);
}

template<typename T, size_t dim> 
mdarray<T,dim>::mdarray(T* Ptr_, const std::array<size_t,dim>& Size_)
{
    this->initialize(Ptr_, Size_);
}

template<typename T, size_t dim> 
void mdarray<T,dim>::TotalSizeAndOffset()
{
    TotalSize = 1;
    for(int idim=0; idim<dim; idim++){
        TotalSize *= Size[idim]; 
    }
    Offset[dim-1] = 1;
    for(int idim=int(dim)-2; idim>=0; idim--){
        Offset[idim] = Size[idim+1]*Offset[idim+1];
    }
}

template<typename T, size_t dim> 
void mdarray<T,dim>::initialize(const std::array<size_t,dim>& Size_)
{
    (*this).~mdarray<T,dim>();
    Size = Size_;
    TotalSizeAndOffset();

    //allocate contiguous memory ot totalsize
    this->Ptr = new T[TotalSize];
    std::fill((*this).begin(), (*this).end(), 0.);

}

template<typename T, size_t dim> 
void mdarray<T,dim>::initialize(T* Ptr_, const std::array<size_t,dim>& Size_)
{
    Ptr = Ptr_; 
    Size = Size_; 
    NotDestruct=true;    //The Ptr_ must be destructed by the other owners. care!
 
    TotalSizeAndOffset();
}

template <typename T, size_t dim>
void mdarray<T,dim>::fill(const T& FillingValue)
{
    std::fill((*this).begin(), (*this).end(), FillingValue);
}

template <typename T, size_t dim>
template <typename... Args>
int mdarray<T,dim>::oneDindex(Args... args) const
{
    assert(dim == sizeof...(args));
    std::array<int ,dim> i = {args...};
    for(int idim=0; idim<dim; idim++){
        assert(i[idim] >=0 && i[idim] < Size[idim]);
    }
    int oneDindex = i[dim-1]*Offset[dim-1];
    for(int idim=int(dim)-2; idim>=0; idim--){
        oneDindex += Offset[idim]*i[idim];               
    }
    return oneDindex;
}


template <typename T, size_t dim>
std::array<int,dim> mdarray<T,dim>::nDindex(const auto& oneDindex) const
{
    std::array<int,dim> nDindex;
    int remainder = oneDindex;
    for(int idim=dim-1; idim>=0; idim--){
        nDindex[idim] = remainder/Offset[idim];
        remainder = remainder%Offset[idim];
    }
}


template <typename T, size_t dim>
template <typename... Args>
T const& mdarray<T,dim>::operator()(Args... args) const 
{
    return Ptr[oneDindex(args...)]; 
}


template <typename T, size_t dim>
template <typename... Args>
inline T& mdarray<T,dim>::operator()(Args... args)
{
    return (const_cast<T&>(static_cast<mdarray<T,dim> const&>(*this)(args...)));
}


template <typename T, size_t dim>
inline T const& mdarray<T,dim>::operator[](const int& oneDindex) const
{
    assert(oneDindex>=0 && oneDindex < TotalSize);
    return Ptr[oneDindex]; 
}


template <typename T, size_t dim>
inline T& mdarray<T,dim>::operator[](const int& oneDindex) 
{
    return (const_cast<T&>(static_cast<mdarray<T,dim> const&>(*this)[oneDindex]));
}


template <typename T, size_t dim>
inline const size_t mdarray<T,dim>::get_Size(const int& index) const
{
    return Size[index];
}

template <typename T, size_t dim>
mdarray<T,dim>::~mdarray()
{
    if(this->Ptr != nullptr && !NotDestruct){
        delete[] Ptr;
    }
    Ptr=nullptr;
}

template<typename T, size_t dim> //requires ( dim>0 && dim<7 )
std::ostream& operator<<(std::ostream& os, const mdarray<T,dim>& mdarray_)
{
    for(auto& number : mdarray_){
        os << number;
        os << std::endl;
    }
    return os;
}
