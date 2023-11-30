#ifndef MDCONTAINERS_HPP
#define MDCONTAINERS_HPP

#include <iostream>
#include <memory>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <memory>

#include <iterator> // For std::forward_iterator_tag
#include <cstddef>  // For std::ptrdiff_t

template<typename T, size_t dim> //requires ( dim>0 && dim<7 )
class mdarray
{
    private:
        //std::unique_ptr<T[]> Ptr =nullptr;
        T* Ptr=nullptr;
        std::array<size_t, dim> Size{0};
        size_t TotalSize=0;
        std::array<size_t, dim> Offset{0};
        bool NotDestruct = false;
        void TotalSizeAndOffset();

    public:
        mdarray(){Ptr = nullptr;};
        mdarray(const mdarray<T,dim>& ToBeCopied);
        mdarray<T,dim>& operator=(const mdarray<T,dim>& ToBeCopied);
        
        mdarray(mdarray<T,dim>&& ToBeMoved);
        mdarray<T,dim>& operator=(mdarray<T,dim>&& ToBeMoved);

        mdarray(const std::array<size_t,dim>& Size_);        
        void initialize(const std::array<size_t,dim>& Size_);

        mdarray(T* Ptr_, const std::array<size_t,dim>& Size_);
        void initialize(T* Ptr_, const std::array<size_t,dim>& Size_);

        struct Iterator
        {
            //iterator tags
            using iterator_category = std::forward_iterator_tag;
            using difference_type   = std::ptrdiff_t;
            using value_type        = T;
            using pointer           = value_type*; 
            using reference         = value_type&;  

            //iterator constructor
            Iterator(pointer ptr) : m_ptr(ptr) {}
        
            //iterator overloadings
            reference operator*() const { return *m_ptr; }
            pointer operator->() { return m_ptr; }
            // Prefix increment
            Iterator& operator++() { m_ptr++; return *this; }  
            // Postfix increment
            Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
            friend bool operator== (const Iterator& a, const Iterator& b) { return a.m_ptr == b.m_ptr; };
            friend bool operator!= (const Iterator& a, const Iterator& b) { return a.m_ptr != b.m_ptr; };     
        
        private:
            pointer m_ptr;
        };

        Iterator begin() const{ return Iterator(Ptr); }
        Iterator end() const{ return Iterator(Ptr+TotalSize); } // TotalSize is out of bounds        
        const auto& data() const {return Ptr;};
        auto& data() {return Ptr;};
        
        template <typename... Args>
        inline size_t oneDindex(Args... args) const;
        
        inline std::array<size_t,dim> nDindex(const auto& oneDindex) const;

        template <typename... Args>
        inline T const& operator()(Args... args) const;

        template <typename... Args>
        inline T& operator()(Args... args);
        
        inline T const& operator[](const int& oneDindex) const;
        inline T& operator[](const int& oneDindex);        


        inline const size_t get_Size(const int& index) const;
        ~mdarray();

        template<typename T_, size_t dim_>
        friend std::ostream& operator<<(std::ostream&, const mdarray<T_,dim_>& mdarray_); 
};

#include "mdContainers.cpp"


#endif