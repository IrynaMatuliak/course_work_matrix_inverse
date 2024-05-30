#pragma once
#include <iostream>
#include <iomanip>
#include <mem.h>
#include <cassert>
#include "matrixinverse.h"

//#define _TEMPLATE_SPECIFICATION_

class SquareMatrix {
public:
    SquareMatrix() = delete;

    explicit SquareMatrix(size_t n) : dimension(n) {
        data = new double[dimension * dimension];
        init();
    }

    SquareMatrix(const SquareMatrix& m) : dimension(m.dimension), data(new double[dimension * dimension]) {
        init();
    }

    ~SquareMatrix() {
        delete[] data;
    }

    size_t size() {
        return dimension;
    }

    SquareMatrix inverse(int method) {
        assert(method >=0 && method <= 1);
        return (_inverse[method]->*ptr)(*this);
        //return _inverse[method]->inverse(*this);  // The same as above
    }

    double& operator()(size_t col, size_t row) {
        return data[row * dimension + col];
    }
    const double& operator()(size_t col, size_t row) const {
        return data[row * dimension + col];
    }

private:
    class RowProxy
    {
    private:
        size_t offset;
        double* dataPtr;
    public:
        RowProxy() = delete;
        explicit RowProxy(double* ptr, size_t offset) : dataPtr(ptr), offset(offset) {}
        double& operator[](size_t j)
        {
            return dataPtr[j * offset];
        }
        const double& operator[](size_t j) const
        {
            return dataPtr[j * offset];
        }
    };

public:
    RowProxy operator[](size_t i) const {
        return RowProxy(data + i, dimension);
    }

    //friend std::ostream& operator<<(std::ostream& out, const SquareMatrix& matrix);

    friend SquareMatrix CLUPMatrixInverse::inverse(const SquareMatrix&);
    friend SquareMatrix CShultsMatrixInverse::inverse(const SquareMatrix&);

private:
    size_t dimension;
    double* data;
    CLUPMatrixInverse _LUPInverse;
    CShultsMatrixInverse _ShultsInverse;
    CBaseMatrixInverse* _inverse[2];
    SquareMatrix(CBaseMatrixInverse::* ptr)(const SquareMatrix&) = &CBaseMatrixInverse::inverse;

    void init(void)
    {
        _inverse[0] = &_LUPInverse;
        _inverse[1] = &_ShultsInverse;

        memset(data, 0, sizeof(double) * dimension * dimension);
        /*for (size_t i = 0; i < dimension; ++i) {
            for (size_t j = 0; j < dimension; ++j) {
                data[i * dimension + j] = static_cast<double>(0.);
            }
        }*/
    }
};



/*std::ostream& operator<<(std::ostream& out, const SquareMatrix& matrix)
{
    for (size_t i = 0; i < matrix.dimension; ++i) {
        for (size_t j = 0; j < matrix.dimension; ++j) {
#ifndef _TEMPLATE_SPECIFICATION_
                out << std::setw(10) << std::fixed << std::setprecision(4) << matrix.data[i * matrix.dimension + j];
#else
            out << std::setw(10) << matrix.data[i * matrix.dimension + j];
#endif
        }
        out << "\n";
    }
    out << std::endl;
    return out;
}*/

#ifdef _TEMPLATE_SPECIFICATION_
template <>
std::ostream& operator<<(std::ostream& out, const SquareMatrix<double>& matrix)
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
#endif  // Template specification
