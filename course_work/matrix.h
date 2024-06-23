#pragma once
#include <iostream>
#include <vector>
#include <optional>
#include "matrixinverse.h"

class SquareMatrix {
public:
    SquareMatrix() = delete;
    explicit SquareMatrix(size_t n);
    SquareMatrix(const std::vector<std::vector<long double>>& v);
    SquareMatrix(const SquareMatrix& m);
    SquareMatrix(SquareMatrix&& m) noexcept;
    SquareMatrix& operator=(const SquareMatrix& m);
    SquareMatrix& operator=(SquareMatrix&& m);
    ~SquareMatrix();

    size_t size() const;
    Results inverse(int method, long double precision = 0.);
    long double& operator()(size_t col, size_t row);
    const long double& operator()(size_t col, size_t row) const;
    SquareMatrix operator*(const SquareMatrix& m) const;
    SquareMatrix operator*(long double val) const;
    SquareMatrix operator+(const SquareMatrix& m) const;
    void transpond();
    void normalize();
    long double det() const;

    static SquareMatrix IdentityMatrix(size_t n);
    void permutate_rows(size_t r1, size_t r2);
    void permutate_columns(size_t c1, size_t c2);
    SquareMatrix LowTriangularMatrix();
    SquareMatrix HighTriangularMatrix();
    bool isSingular();

private:
    class RowProxy
    {
    private:
        size_t offset;
        long double* dataPtr;
    public:
        RowProxy() = delete;
        explicit RowProxy(long double* ptr, size_t offset) : dataPtr(ptr), offset(offset) {}
        long double& operator[](size_t j)
        {
            return dataPtr[j * offset];
        }
        const long double& operator[](size_t j) const
        {
            return dataPtr[j * offset];
        }
    };

public:
    RowProxy operator[](size_t i) const {
        return RowProxy(data + i, dimension);
    }

    friend std::ostream& operator<<(std::ostream& out, const SquareMatrix& matrix);
    friend Results CLUPMatrixInverse::inverse(const SquareMatrix&, double precision);
    friend Results CShultsMatrixInverse::inverse(const SquareMatrix&, double precision);

private:
    size_t dimension;
    long double* data;
    CLUPMatrixInverse _LUPInverse;
    CShultsMatrixInverse _ShultsInverse;
    CBaseMatrixInverse* _inverse[2];
    Results(CBaseMatrixInverse::* ptr)(const SquareMatrix&, double) = &CBaseMatrixInverse::inverse;

    void init(void);
    void init(const SquareMatrix& m);
    std::optional<size_t> find_zero() const;
    SquareMatrix prepare_for_Gaus(long double& sign) const;
};

struct Results {
    SquareMatrix A;
    SquareMatrix LU;
    SquareMatrix X;
    SquareMatrix P;
    long int complexity;
    long double precision;
    int method;
    long double determinant;

    Results();
};

