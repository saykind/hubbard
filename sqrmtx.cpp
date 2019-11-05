#include "sqrmtx.h"
#include <iostream>
#include <iomanip>

/* SqrMtx class */
// Initialization
template <typename data_t> 
SqrMtx<data_t>::SqrMtx(int N_, const data_t& d) {
    if (N_ < 1)
        throw std::invalid_argument("non-positive matrix size.");
    N = N_;
    x = new data_t*[N];
    for (int i = 0; i<N; i++) {
        x[i] = new data_t[N];
        for (int j=0; j<N; j++)
            x[i][j] = 0;
        x[i][i] = d;
    }
}

template <typename data_t> 
SqrMtx<data_t>::SqrMtx(const SqrMtx<data_t>& m) {
    N = m.N;
    x = new data_t*[N];
    for (int i = 0; i<N; i++) {
        x[i] = new data_t[N];
        for (int j=0; j<N; j++)
            x[i][j] = m.x[i][j];
    }
}

template <typename data_t> 
SqrMtx<data_t>::~SqrMtx() {
    if (!x) {return;}
    for (int i = 0; i<N; i++)
        delete[] x[i];
    delete[] x;
    x = nullptr;
}

// Methods
template <typename data_t>
SqrMtx<data_t> SqrMtx<data_t>::adj(int i0, int j0) const {
    if (!((0 <= i0 < N) && (0 <= j0 < N)))
        throw std::invalid_argument("matrix index is out of range.");
    if (N < 2)
        throw std::invalid_argument("adjugate matrix cannot be calculated.");
    SqrMtx<data_t> a(N-1);
    for (int i=0; i+1<N; i++) {
        for (int j=0; j+1<N; j++)
            a[i][j] = x[i<i0 ? i : i+1][j<j0 ? j : j+1];
    }
    return a;
}

template <typename data_t>
data_t SqrMtx<data_t>::det() const { // TODO implement faster algorithm
    if (N < 2)
        return x[0][0];
    data_t d = 0;
    for (int i=0; i<N; i++)
        d += (i%2 ? -1 : 1)*x[i][0]*adj(i,0).det();
    return d;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::swp_rws(int i1, int i2) {
    for (int j=0; j<N; j++) {
        data_t swp = x[i1][j];
        x[i1][j] = x[i2][j];
        x[i2][j] = swp;
    }
    return *this;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::inv() {
    SqrMtx<data_t> r(N,1);
    for (int j=0; j+1<N; j++) {
        if (!x[j][j])
            for (int i=j+1; i<N; i++) {
                if (x[i][j]) {
                    swp_rws(j,i);
                    r.swp_rws(j,i);
                    break;
                }
                if (i+1==N)
                    throw std::invalid_argument("non-invertable matrix");
            }
        for (int i=j+1; i<N; i++) {
            data_t f = x[i][j]/x[j][j];
            for (int k=0; k<N; k++) {
                x[i][k] -= x[j][k]*f;
                r[i][k] -= r[j][k]*f;
            }
        }
    }
    for (int j=N-1; j>0; j--) {
        for (int i=j-1; i>=0; i--) {
            data_t f = x[i][j]/x[j][j];
            for (int k=N-1; k>=0; k--) {
                x[i][k] -= x[j][k]*f;
                r[i][k] -= r[j][k]*f;
            }
        }
    }
    for (int i=0; i<N; i++) {
        data_t f = x[i][i];
        for (int j=0; j<N; j++)
            x[i][j] = r[i][j]/f;
    }
    return *this;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::chp(data_t lmt) {
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            x[i][j] = abs(x[i][j]) > lmt ? x[i][j] : 0;
    return *this;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::rnd(int max) {
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            x[i][j] = rand()%max;
    return *this;
}

// Overloaded member-operators
template <typename data_t>
data_t* SqrMtx<data_t>::operator[](int i) {return x[i];}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator=(const SqrMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            x[i][j] = m.x[i][j];
    return *this;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator=(const DglMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++)
            x[i][j] = 0;
        x[i][i] = m.x[i];
    }
    return *this;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator=(const data_t& d) {
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++)
            x[i][j] = 0;
        x[i][i] = d;
    }
    return *this;
}

// Arithmetic operators
template <typename data_t>
SqrMtx<data_t> SqrMtx<data_t>::operator+(const data_t& d) const {
    SqrMtx<data_t> r(*this);
    for (int i=0; i<N; i++)
        r.x[i][i] += d;
    return r;
}

template <typename data_t>
SqrMtx<data_t> SqrMtx<data_t>::operator-(const data_t& d) const {
    SqrMtx<data_t> r(*this);
    for (int i=0; i<N; i++)
        r.x[i][i] -= d;
    return r;
}

template <typename data_t>
SqrMtx<data_t> SqrMtx<data_t>::operator*(const data_t& d) const {
    SqrMtx<data_t> r(*this);
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            r.x[i][j] *= d;
    return r;
}

template <typename data_t>
SqrMtx<data_t> SqrMtx<data_t>::operator/(const data_t& d) const {
    SqrMtx<data_t> r(*this);
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            r.x[i][j] /= d;
    return r;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator+=(const data_t& d) {
    for (int i=0; i<N; i++)
        x[i][i] += d;
    return *this;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator-=(const data_t& d) {
    for (int i=0; i<N; i++)
        x[i][i] -= d;
    return *this;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator*=(const data_t& d) {
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            x[i][j] *= d;
    return *this;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator/=(const data_t& d) {
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            x[i][j] /= d;
    return *this;
}

// Arithemetic with SqrMtx
template <typename data_t>
SqrMtx<data_t> SqrMtx<data_t>::operator+(const SqrMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    SqrMtx<data_t> r(N);
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            r.x[i][j] = x[i][j] + m.x[i][j];
    return r;
}

template <typename data_t>
SqrMtx<data_t> SqrMtx<data_t>::operator-(const SqrMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    SqrMtx<data_t> r(N);
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            r.x[i][j] = x[i][j] - m.x[i][j];
    return r;
}

template <typename data_t>
SqrMtx<data_t> SqrMtx<data_t>::operator*(const SqrMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    SqrMtx<data_t> r(N);
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++) {
            data_t sum = 0;
            for (int l=0; l<N; l++)
                sum += x[i][l]*m.x[l][j];
            r[i][j] = sum;
        }
    return r;
}

template <typename data_t>
SqrMtx<data_t> SqrMtx<data_t>::operator/(const SqrMtx<data_t>& m) const {
    SqrMtx<data_t> r(m);
    return *this*r.inv();
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator+=(const SqrMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            x[i][j] += m.x[i][j];
    return *this;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator-=(const SqrMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            x[i][j] -= m.x[i][j];
    return *this;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator*=(const SqrMtx<data_t>& m) {
    return *this = *this*m;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator/=(const SqrMtx<data_t>& m) {
    SqrMtx<data_t> r(m);
    return *this = *this*r.inv();
}

// Arithmetic operators with DglMtx
template <typename data_t>
SqrMtx<data_t> SqrMtx<data_t>::operator+(const DglMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    SqrMtx<data_t> r(N);
    for (int i=0; i<N; i++)
        r.x[i][i] = x[i][i] + m.x[i];
    return r;
}

template <typename data_t>
SqrMtx<data_t> SqrMtx<data_t>::operator-(const DglMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    SqrMtx<data_t> r(N);
    for (int i=0; i<N; i++)
        r.x[i][i] = x[i][i] - m.x[i];
    return r;
}

template <typename data_t>
SqrMtx<data_t> SqrMtx<data_t>::operator*(const DglMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    SqrMtx<data_t> r(N);
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            r[i][j] = x[i][j]*m.x[j];
    return r;
}

template <typename data_t>
SqrMtx<data_t> SqrMtx<data_t>::operator/(const DglMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    SqrMtx<data_t> r(N);
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            r[i][j] = x[i][j]/m.x[j];
    return r;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator+=(const DglMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++)
        x[i][i] += m.x[i];
    return *this;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator-=(const DglMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++)
        x[i][i] -= m.x[i];
    return *this;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator*=(const DglMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            x[i][j] = x[i][j]*m.x[j];
    return *this;
}

template <typename data_t>
SqrMtx<data_t>& SqrMtx<data_t>::operator/=(const DglMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            x[i][j] = x[i][j]/m.x[j];
    return *this;
}


// Overloaded friend-operators
template <typename data_t>
std::ostream& operator<<(std::ostream& stream, const SqrMtx<data_t>& m) {
    stream << '[';
    for (int i=0; i<m.N; i++) {
        if (i) {stream << ' ';}
        stream << "[  ";
        for (int j=0; j<m.N; j++)
           stream << std::setw(6) << m.x[i][j] << "  ";
        stream << ']';
        if ((i+1)<m.N) {stream << std::endl;}
    }
    return stream << ']' << std::endl;
}

template <typename data_t>
SqrMtx<data_t> operator+(const data_t& s, const SqrMtx<data_t>& m) {
    SqrMtx<data_t> r(m);
    for (int i=0; i<r.N; i++)
        r[i][i] += s;
    return r;
}

template <typename data_t>
SqrMtx<data_t> operator-(const data_t& s, const SqrMtx<data_t>& m) {
    SqrMtx<data_t> r(m);
    for (int i=0; i<r.N; i++){
        r[i][i] -= s;
        for (int j=0; j<r.N; j++)
            r[i][j] *= -1;
    }
    return r;
}

template <typename data_t>
SqrMtx<data_t> operator*(const data_t& s, const SqrMtx<data_t>& m) {
    SqrMtx<data_t> r(m);
    for (int i=0; i<r.N; i++)
        for (int j=0; j<r.N; j++)
            r[i][j] *= s;
    return r;
}

template <typename data_t>
SqrMtx<data_t> operator/(const data_t& s, const SqrMtx<data_t>& m) {
    SqrMtx<data_t> r(m);
    return s*r.inv();
}











/* DglMtx class */
// Initialization
template <typename data_t> 
DglMtx<data_t>::DglMtx(int N_, const data_t& d) {
    if (N_ < 1)
        throw std::invalid_argument("non-positive matrix size.");
    N = N_;
    x = new data_t[N];
    for (int i = 0; i<N; i++)
        x[i] = d;
}

template <typename data_t> 
DglMtx<data_t>::DglMtx(const DglMtx<data_t>& m) {
    N = m.N;
    x = new data_t[N];
    for (int i = 0; i<N; i++)
        x[i] = m.x[i];
}

template <typename data_t> 
DglMtx<data_t>::~DglMtx() {
    if (!x) {return;}
    delete[] x;
    x = nullptr;
}

// Overloaded member-operators
template <typename data_t>
data_t& DglMtx<data_t>::operator[](int i) {return x[i];}

template <typename data_t>
DglMtx<data_t>& DglMtx<data_t>::operator=(const DglMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++)
        x[i] = m.x[i];
    return *this;
}

template <typename data_t>
DglMtx<data_t>& DglMtx<data_t>::operator=(const SqrMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++)
        x[i] = m.x[i][i];
    return *this;
}

template <typename data_t>
DglMtx<data_t>& DglMtx<data_t>::operator=(const data_t& d) {
    for (int i=0; i<N; i++) {
        x[i] = d;
    }
    return *this;
}

// Arithmetic operators with scalars
template <typename data_t>
DglMtx<data_t> DglMtx<data_t>::operator+(const data_t& d) const {
    DglMtx<data_t> r(*this);
    for (int i=0; i<N; i++)
        r.x[i] += d;
    return r;
}

template <typename data_t>
DglMtx<data_t> DglMtx<data_t>::operator-(const data_t& d) const {
    DglMtx<data_t> r(*this);
    for (int i=0; i<N; i++)
        r.x[i] -= d;
    return r;
}

template <typename data_t>
DglMtx<data_t> DglMtx<data_t>::operator*(const data_t& d) const {
    DglMtx<data_t> r(*this);
    for (int i=0; i<N; i++)
        r.x[i] *= d;
    return r;
}

template <typename data_t>
DglMtx<data_t> DglMtx<data_t>::operator/(const data_t& d) const {
    DglMtx<data_t> r(*this);
    for (int i=0; i<N; i++)
        r.x[i] /= d;
    return r;
}

template <typename data_t>
DglMtx<data_t>& DglMtx<data_t>::operator+=(const data_t& d) {
    for (int i=0; i<N; i++)
        x[i] += d;
    return *this;
}

template <typename data_t>
DglMtx<data_t>& DglMtx<data_t>::operator-=(const data_t& d) {
    for (int i=0; i<N; i++)
        x[i] -= d;
    return *this;
}

template <typename data_t>
DglMtx<data_t>& DglMtx<data_t>::operator*=(const data_t& d) {
    for (int i=0; i<N; i++)
        x[i] *= d;
    return *this;
}

template <typename data_t>
DglMtx<data_t>& DglMtx<data_t>::operator/=(const data_t& d) {
    for (int i=0; i<N; i++)
        x[i] /= d;
    return *this;
}

// Arithmetic operators with DglMtx
template <typename data_t>
DglMtx<data_t> DglMtx<data_t>::operator+(const DglMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    DglMtx<data_t> r(N);
    for (int i=0; i<N; i++)
        r.x[i] = x[i] + m.x[i];
    return r;
}

template <typename data_t>
DglMtx<data_t> DglMtx<data_t>::operator-(const DglMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    DglMtx<data_t> r(N);
    for (int i=0; i<N; i++)
        r.x[i] = x[i] - m.x[i];
    return r;
}

template <typename data_t>
DglMtx<data_t> DglMtx<data_t>::operator*(const DglMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    DglMtx<data_t> r(N);
    for (int i=0; i<N; i++)
        r[i] = x[i]*m.x[i];
    return r;
}

template <typename data_t>
DglMtx<data_t> DglMtx<data_t>::operator/(const DglMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    DglMtx<data_t> r(N);
    for (int i=0; i<N; i++)
        r[i] = x[i]/m.x[i];
    return r;
}

template <typename data_t>
DglMtx<data_t>& DglMtx<data_t>::operator+=(const DglMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++)
        x[i] += m.x[i];
    return *this;
}

template <typename data_t>
DglMtx<data_t>& DglMtx<data_t>::operator-=(const DglMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++)
        x[i] -= m.x[i];
    return *this;
}

template <typename data_t>
DglMtx<data_t>& DglMtx<data_t>::operator*=(const DglMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++)
        x[i] *= m.x[i];
    return *this;
}

template <typename data_t>
DglMtx<data_t>& DglMtx<data_t>::operator/=(const DglMtx<data_t>& m) {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    for (int i=0; i<N; i++)
        x[i] /= m.x[i];
    return *this;
}

// Arithmetic operators with SqrMtx
template <typename data_t>
SqrMtx<data_t> DglMtx<data_t>::operator+(const SqrMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    SqrMtx<data_t> r(m);
    for (int i=0; i<N; i++)
        r.x[i][i] += x[i];
    return r;
}

template <typename data_t>
SqrMtx<data_t> DglMtx<data_t>::operator-(const SqrMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    SqrMtx<data_t> r(m);
    for (int i=0; i<N; i++) {
        r.x[i][i] -= x[i];
        for (int j=0; j<N; i++)
            r.x[i][j] *= -1;
    }
    return r;
}

template <typename data_t>
SqrMtx<data_t> DglMtx<data_t>::operator*(const SqrMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    SqrMtx<data_t> r(m);
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            r[i][j] = x[i]*m.x[i][j];
    return r;
}

template <typename data_t>
SqrMtx<data_t> DglMtx<data_t>::operator/(const SqrMtx<data_t>& m) const {
    if (N != m.N)
        throw std::invalid_argument("mismatched matrix size.");
    SqrMtx<data_t> r(m);
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            r[i][j] = x[i]/m.x[i][j];
    return r;
}

// Overloaded friend-operators
template <typename data_t>
std::ostream& operator<<(std::ostream& stream, const DglMtx<data_t>& m) {
    stream << "[[";
    for (int i=0; i<m.N; i++) {
        stream << std::setw(6) << m.x[i] << "  ";
    }
    return stream << "]]" << std::endl;
}

template <typename data_t>
DglMtx<data_t> operator+(const data_t& s, const DglMtx<data_t>& m) {
    DglMtx<data_t> r(m);
    for (int i=0; i<r.N; i++)
        r[i] += s;
    return r;
}

template <typename data_t>
DglMtx<data_t> operator-(const data_t& s, const DglMtx<data_t>& m) {
    DglMtx<data_t> r(m);
    for (int i=0; i<r.N; i++)
        r[i] = s-r[i];
    return r;
}

template <typename data_t>
DglMtx<data_t> operator*(const data_t& s, const DglMtx<data_t>& m) {
    DglMtx<data_t> r(m);
    for (int i=0; i<r.N; i++)
        r[i] *= s;
    return r;
}

template <typename data_t>
DglMtx<data_t> operator/(const data_t& s, const DglMtx<data_t>& m) {
    DglMtx<data_t> r(m);
    for (int i=0; i<r.N; i++)
        r[i] *= s/r[i];
    return r;
}









/* Template implementations */
// SqrMtx class
template class SqrMtx<int>;
template std::ostream& operator<<<int>(std::ostream& stream, const SqrMtx<int>& m);
template SqrMtx<int> operator+<int>(const int& s, const SqrMtx<int>& m);
template SqrMtx<int> operator-<int>(const int& s, const SqrMtx<int>& m);
template SqrMtx<int> operator*<int>(const int& s, const SqrMtx<int>& m);
template SqrMtx<int> operator/<int>(const int& s, const SqrMtx<int>& m);

template class SqrMtx<double>;
template std::ostream& operator<<<double>(std::ostream& stream, const SqrMtx<double>& m);
template SqrMtx<double> operator+<double>(const double& s, const SqrMtx<double>& m);
template SqrMtx<double> operator-<double>(const double& s, const SqrMtx<double>& m);
template SqrMtx<double> operator*<double>(const double& s, const SqrMtx<double>& m);
template SqrMtx<double> operator/<double>(const double& s, const SqrMtx<double>& m);

/* DglMtx */
template class DglMtx<int>;
template std::ostream& operator<<<int>(std::ostream& stream, const DglMtx<int>& m);
template DglMtx<int> operator+<int>(const int& s, const DglMtx<int>& m);
template DglMtx<int> operator-<int>(const int& s, const DglMtx<int>& m);
template DglMtx<int> operator*<int>(const int& s, const DglMtx<int>& m);
template DglMtx<int> operator/<int>(const int& s, const DglMtx<int>& m);

template class DglMtx<double>;
template std::ostream& operator<<<double>(std::ostream& stream, const DglMtx<double>& m);
template DglMtx<double> operator+<double>(const double& s, const DglMtx<double>& m);
template DglMtx<double> operator-<double>(const double& s, const DglMtx<double>& m);
template DglMtx<double> operator*<double>(const double& s, const DglMtx<double>& m);
template DglMtx<double> operator/<double>(const double& s, const DglMtx<double>& m);

