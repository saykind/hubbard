#pragma once
#include <iostream>
#include "defaults.h"

template <typename data_t> 
class SqrMtx;

template <typename data_t> 
class DglMtx;

typedef SqrMtx<double>  mtx;
typedef SqrMtx<int>     mxt;

typedef DglMtx<double>  dmx;
typedef DglMtx<int>     dtx;

template <typename data_t>
class SqrMtx {
        template <typename type_t>
        friend class DglMtx;
    private:
        int N = 0;
        data_t** x = nullptr;
    public:
                SqrMtx(int N_, const data_t& d = 0);
                SqrMtx(const SqrMtx&);
                ~SqrMtx();
    public:
        SqrMtx  adj(int, int) const;
        data_t  det() const;
        SqrMtx& inv();
        SqrMtx& swp_rws(int, int);
        SqrMtx& chp(data_t lmt = 1e-14);
        SqrMtx& rnd(int max = RAND_MAX);
    public:
        data_t* operator[](int);
        SqrMtx& operator=(const SqrMtx&);
        SqrMtx& operator=(const DglMtx<data_t>&);
        SqrMtx& operator=(const data_t&);
    public:
        SqrMtx  operator+(const data_t&) const;
        SqrMtx  operator-(const data_t&) const;
        SqrMtx  operator*(const data_t&) const;
        SqrMtx  operator/(const data_t&) const;
        SqrMtx& operator+=(const data_t&);
        SqrMtx& operator-=(const data_t&);
        SqrMtx& operator*=(const data_t&);
        SqrMtx& operator/=(const data_t&);
    public:
        SqrMtx  operator+(const SqrMtx&) const;
        SqrMtx  operator-(const SqrMtx&) const;
        SqrMtx  operator*(const SqrMtx&) const;
        SqrMtx  operator/(const SqrMtx&) const;
        SqrMtx& operator+=(const SqrMtx&);
        SqrMtx& operator-=(const SqrMtx&);
        SqrMtx& operator*=(const SqrMtx&);
        SqrMtx& operator/=(const SqrMtx&);
    public:
        SqrMtx  operator+(const DglMtx<data_t>&) const;
        SqrMtx  operator-(const DglMtx<data_t>&) const;
        SqrMtx  operator*(const DglMtx<data_t>&) const;
        SqrMtx  operator/(const DglMtx<data_t>&) const;
        SqrMtx& operator+=(const DglMtx<data_t>&);
        SqrMtx& operator-=(const DglMtx<data_t>&);
        SqrMtx& operator*=(const DglMtx<data_t>&);
        SqrMtx& operator/=(const DglMtx<data_t>&);
    public:
        template <typename type_t>
        friend std::ostream& operator<<(std::ostream&, const SqrMtx<type_t>&);
        template <typename type_t>
        friend SqrMtx<type_t> operator+(const type_t&, const SqrMtx<type_t>&);
        template <typename type_t>
        friend SqrMtx<type_t> operator-(const type_t&, const SqrMtx<type_t>&);
        template <typename type_t>
        friend SqrMtx<type_t> operator*(const type_t&, const SqrMtx<type_t>&);
        template <typename type_t>
        friend SqrMtx<type_t> operator/(const type_t&, const SqrMtx<type_t>&);
};

template <typename data_t>
class DglMtx {
        template <typename type_t>
        friend class SqrMtx;
    private:
        int N = 0;
        data_t* x;
    public:
                DglMtx(int N, const data_t& d = 0);
                DglMtx(const DglMtx& m);
                ~DglMtx();
    public:
        DglMtx  adj(int, int) const;
        data_t  det() const;
        DglMtx  inv() const;
        DglMtx  chp(data_t lmt = 1e-14);
        DglMtx& rnd(int max = RAND_MAX);
    public:
        data_t& operator[](int);
        DglMtx& operator=(const DglMtx&);
        DglMtx& operator=(const SqrMtx<data_t>&);
        DglMtx& operator=(const data_t&);
    public:
        DglMtx  operator+(const data_t&) const;
        DglMtx  operator-(const data_t&) const;
        DglMtx  operator*(const data_t&) const;
        DglMtx  operator/(const data_t&) const;
        DglMtx& operator+=(const data_t&);
        DglMtx& operator-=(const data_t&);
        DglMtx& operator*=(const data_t&);
        DglMtx& operator/=(const data_t&);
    public:
        DglMtx  operator+(const DglMtx&) const;
        DglMtx  operator-(const DglMtx&) const;
        DglMtx  operator*(const DglMtx&) const;
        DglMtx  operator/(const DglMtx&) const;
        DglMtx& operator+=(const DglMtx&);
        DglMtx& operator-=(const DglMtx&);
        DglMtx& operator*=(const DglMtx&);
        DglMtx& operator/=(const DglMtx&);
    public:
        SqrMtx<data_t>  operator+(const SqrMtx<data_t>&) const;
        SqrMtx<data_t>  operator-(const SqrMtx<data_t>&) const;
        SqrMtx<data_t>  operator*(const SqrMtx<data_t>&) const;
        SqrMtx<data_t>  operator/(const SqrMtx<data_t>&) const;
    public:
        template <typename type_t>
        friend std::ostream& operator<<(std::ostream&, const DglMtx<type_t>&);
        template <typename type_t>
        friend DglMtx<type_t> operator+(const type_t&, const DglMtx<type_t>&);
        template <typename type_t>
        friend DglMtx<type_t> operator-(const type_t&, const DglMtx<type_t>&);
        template <typename type_t>
        friend DglMtx<type_t> operator*(const type_t&, const DglMtx<type_t>&);
        template <typename type_t>
        friend DglMtx<type_t> operator/(const type_t&, const DglMtx<type_t>&);
};

