#pragma once
#include <iostream>
#include <iomanip>
#include <math.h>
//#include <mem.h>
#include <cassert>
#include "matrix.h"

SquareMatrix::SquareMatrix(size_t n) : dimension(n), data(new double[dimension * dimension])
{
    init();
}

SquareMatrix::SquareMatrix(const SquareMatrix& m)
            : dimension(m.dimension), data(new double[dimension * dimension]
           /*, _LUPInverse(m._LUPInverse), _ShultsInversem(m._ShultsInverse)*/
           )
{
    init(m);
}

SquareMatrix::SquareMatrix(const std::vector<std::vector<double>>& v)
{
    size_t s = v.size();
    bool f = true;
    for (size_t i = 0; i < s; ++i) {
        if (v[i].size() != s) {
            f = false;
            break;
        }
    }
    if (f) {
        dimension = s;
        data = new double[dimension * dimension];
        size_t c = 0;
        for (size_t j = 0; j < dimension; ++j) {
            for (size_t i = 0; i < dimension; ++i) {
                data[c++] = v[j][i];
            }
        }
        _inverse[0] = &_LUPInverse;
        _inverse[1] = &_ShultsInverse;
    } else {
        dimension = 1;
        data = new double[dimension * dimension];
        init();
    }
}

SquareMatrix::SquareMatrix(SquareMatrix&& m) noexcept
{
    dimension = m.dimension;
    data = new double[1];
    std::swap(m.data, data);
    _inverse[0] = &_LUPInverse;
    _inverse[1] = &_ShultsInverse;
}

SquareMatrix& SquareMatrix::operator=(const SquareMatrix& m)
{
    dimension = m.dimension;
    data = new double[dimension * dimension];
    init(m);
    return *this;
}

SquareMatrix& SquareMatrix::operator=(SquareMatrix&& m)
{
    SquareMatrix new_m = std::move(m);
    std::swap(new_m.dimension, dimension);
    std::swap(new_m.data, data);
    _inverse[0] = &_LUPInverse;
    _inverse[1] = &_ShultsInverse;
    return *this;
}

SquareMatrix::~SquareMatrix()
{
    delete[] data;
}

size_t SquareMatrix::size() const
{
    return dimension;
}

Results SquareMatrix::inverse(int method, double precision)
{
    assert(method >=0 && method <= 1);
    Results res = (_inverse[method]->*ptr)(*this, precision);
    res.A = *this;
    return res;
    //return _inverse[method]->inverse(*this);  // The same as above
}

double& SquareMatrix::operator()(size_t col, size_t row)
{
    return data[row * dimension + col];
}

const double& SquareMatrix::operator()(size_t col, size_t row) const
{
    return data[row * dimension + col];
}

void SquareMatrix::init(void)
{
    _inverse[0] = &_LUPInverse;
    _inverse[1] = &_ShultsInverse;

    //memset(data, 0, sizeof(double) * dimension * dimension);
    for (size_t i = 0; i < dimension; ++i) {
        for (size_t j = 0; j < dimension; ++j) {
            data[i * dimension + j] = static_cast<double>(0.);
        }
    }
}


void SquareMatrix::init(const SquareMatrix& m)
{
    _inverse[0] = &_LUPInverse;
    _inverse[1] = &_ShultsInverse;

    //memcpy(data, m.data, sizeof(double) * dimension * dimension);
    for (size_t i = 0; i < dimension; ++i) {
        for (size_t j = 0; j < dimension; ++j) {
            data[i * dimension + j] = m.data[i * dimension + j];
        }
    }
}

SquareMatrix SquareMatrix::operator*(const SquareMatrix& m) const
{
    SquareMatrix res{dimension};
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            res[j][i] = 0;
            for (int k = 0; k < dimension; k++) {
                res[j][i] += (*this)[k][i] * m[j][k];
            }
        }
    }
    return res;
}

SquareMatrix SquareMatrix::operator*(double val) const
{
    SquareMatrix res = *this;
    for (size_t i = 0; i < dimension * dimension; ++i) {
        res.data[i] *= val;
    }
    return res;
}

SquareMatrix SquareMatrix::operator+(const SquareMatrix& m) const
{
    SquareMatrix res = *this;
    for (int i = 0; i < dimension * dimension; i++) {
        res.data[i] = res.data[i] + m.data[i];
    }
    return res;
}

void SquareMatrix::normalize()
{
    double N1 = 0;
    double Ninf = 0;
    for(size_t row = 0; row < dimension; row++) {
        double colsum = 0;
        double rowsum = 0;
        for(size_t col = 0; col < dimension; col++){
            rowsum += fabs(data[row * dimension + col]);
            colsum += fabs(data[col * dimension + row]);
        }
        N1 = std::max(colsum, N1);
        Ninf = std::max(rowsum, Ninf);
    }

    double norm = 1. / (N1 * Ninf);
    for (size_t i = 0; i < dimension * dimension; ++i) {
        data[i] *= norm;
    }
}

void SquareMatrix::transpond()
{
    for (size_t i = 0; i < dimension; ++i) {
        for (size_t j = i; j < dimension; ++j) {
            std::swap((*this)[i][j], (*this)[j][i]);
        }
    }
}

double SquareMatrix::det() const
{
    SquareMatrix B = *this;

    for (size_t step = 0; step < dimension - 1; step++) {
        for(size_t row = step + 1; row < dimension; row++) {
            double coeff = -B[row][step] / B[step][step];
            for(int col = step; col < dimension; col++) {
                B[row][col] += B[step][col] * coeff;
            }
        }
    }
    double d = 1.;
    for (size_t i = 0; i < dimension; i++) {
        d *= B[i][i];
    }
    return d;
}

SquareMatrix SquareMatrix::IdentityMatrix(size_t n)
{
    SquareMatrix res{n};
    for (int i = 0; i < n; ++i) {
        res[i][i] = 1.;
    }
    return res;
}

void SquareMatrix::permutate(size_t r1, size_t r2)
{
    if (r1 >= dimension || r2 >= dimension) {
        assert(false);
        return;
    }
    for (int i = 0; i < dimension; ++i) {
        double d = data[i + r1 * dimension];
        data[i + r1 * dimension] = data[i + r2 * dimension];
        data[i + r2 * dimension] = d;
    }
}

SquareMatrix SquareMatrix::LowTriangularMatrix()
{
    SquareMatrix res = *this;
    for (int i = 0; i < dimension; ++i) {
        for(int j = i; j < dimension; ++j) {
            if (i != j) {
                res[j][i] = 0.;
            } else {
                res[j][i] = 1.;
            }
        }
    }
    return res;
}

SquareMatrix SquareMatrix::HighTriangularMatrix()
{
    SquareMatrix res = *this;
    for (int i = 0; i < dimension; ++i) {
        for(int j = 0; j < i; ++j) {
            res[j][i] = 0.;
        }
    }
    return res;
}

bool SquareMatrix::isSingular()
{
    return std::fabs(det()) < std::numeric_limits<double>::epsilon();
}

std::ostream& operator<<(std::ostream& out, const SquareMatrix& matrix)
{
    for (size_t i = 0; i < matrix.dimension; ++i) {
        for (size_t j = 0; j < matrix.dimension; ++j) {
            out << std::setw(10) << std::fixed << std::setprecision(4) << matrix.data[i * matrix.dimension + j];
        }
        out << "\n";
    }
    out << std::endl;
    return out;
}

Results::Results() : A(SquareMatrix(1)), LU(SquareMatrix(1)), X(SquareMatrix(1)), P(SquareMatrix(1))
{
}
/*
Results::Results(const Results& r) : LU(r.LU), P(SquareMatrix(r.P))
{

}

Results::Results(Results&& r) : LU(r.LU), P(SquareMatrix(r.P))
{

}

Results& Results::operator=(const Results& r)
{
    LU = r.LU;
    P = r.P;
}

Results& Results::operator=(Results&& r)
{
    LU = r.LU;
    P = r.P;
}
*/
