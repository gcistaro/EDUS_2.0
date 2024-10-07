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

#include "MultiIndex/MultiIndex.hpp"

template<typename T, size_t dim> //requires ( dim>0 && dim<7 )
class mdarray
{
    private:
        //std::unique_ptr<T[]> Ptr =nullptr;
        T* Ptr=nullptr;
        std::array<int, dim> Size{0};
        int TotalSize=0;
        std::array<int, dim> Offset{0};
        MultiIndex<dim> multindex;
        bool NotDestruct = false;
        void TotalSizeAndOffset();

    public:
        mdarray(){Ptr = nullptr;};
        mdarray(const mdarray<T,dim>& ToBeCopied);
        mdarray<T,dim>& operator=(const mdarray<T,dim>& ToBeCopied);
        
        mdarray(mdarray<T,dim>&& ToBeMoved);
        mdarray<T,dim>& operator=(mdarray<T,dim>&& ToBeMoved);

        mdarray(const std::array<int,dim>& Size_);        
        void initialize(const std::array<int,dim>& Size_);

        mdarray(T* Ptr_, const std::array<int,dim>& Size_);
        void initialize(T* Ptr_, const std::array<int,dim>& Size_);
        
        void fill(const T& FillingValue);

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
        
            auto data() { return m_ptr;}
            //iterator overloadings
            reference operator*() const { return *m_ptr; };
            pointer operator->() { return m_ptr; };
            // Prefix increment
            Iterator& operator++() { m_ptr++; return *this; };  
            // Postfix increment
            Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; };
            Iterator operator+=(int rhs) { pointer tmp = this->data(); tmp += rhs; return Iterator(tmp);}
            friend bool operator== (const Iterator& a, const Iterator& b) { return a.m_ptr == b.m_ptr; };
            friend bool operator!= (const Iterator& a, const Iterator& b) { return a.m_ptr != b.m_ptr; };     
            //difference
            int operator-(const Iterator& it2) {return this->m_ptr - it2.m_ptr;}        
            Iterator operator+(const Iterator& it2) {return Iterator(this->m_ptr + it2.m_ptr);}        
            Iterator operator+(const int& i) {return Iterator(this->m_ptr + i);}        
        private:
            pointer m_ptr;
        };

        Iterator begin() const{ return Iterator(Ptr); }
        Iterator end() const{ return Iterator(Ptr+TotalSize); } // TotalSize is out of bounds        
        const auto& data() const {return Ptr;};
        auto& data() {return Ptr;};
        
        template <typename... Args>
        inline int oneDindex(Args... args) const;
        
        inline std::vector<int> nDindex(const int& oneDindex) const;

        template <typename... Args>
        inline T const& operator()(Args... args) const;

        template <typename... Args>
        inline T& operator()(Args... args);
        
        inline T const& operator[](const int& oneDindex) const;
        inline T& operator[](const int& oneDindex);        


        inline const int get_Size(const int& index) const;
        inline auto get_Size() const {return Size;};
        inline auto get_TotalSize() const {return TotalSize;};
        ~mdarray();

        template<typename T_, size_t dim_>
        friend std::ostream& operator<<(std::ostream&, const mdarray<T_,dim_>& mdarray_); 
};

#include "mdContainers_definitions.hpp"


#endif