#pragma once
#include <iostream>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <optional>
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

/*double SquareMatrix::det() const
{
    SquareMatrix B = *this;

    for (size_t step = 0; step < dimension - 1; ++step) {
        for(size_t row = step + 1; row < dimension; ++row) {
            double coeff = -B[row][step] / B[step][step];
            for(int col = step; col < dimension; ++col) {
                B[row][col] += B[step][col] * coeff;
            }
        }
    }
    double d = 1.;
    for (size_t i = 0; i < dimension; ++i) {
        d *= B[i][i];
    }
    return d;
}*/

/*void SquareMatrix::permutateRows(size_t r1, size_t r2) {
    if (r1 >= dimension || r2 >= dimension) {
        assert(false);
        return;
    }
    for (size_t i = 0; i < dimension; ++i) {
        std::swap((*this)(r1, i), (*this)(r2, i));
    }
}

void SquareMatrix::permutateCols(size_t c1, size_t c2) {
    if (c1 >= dimension || c2 >= dimension) {
        assert(false);
        return;
    }
    for (size_t i = 0; i < dimension; ++i) {
        std::swap((*this)(i, c1), (*this)(i, c2));
    }
}

double SquareMatrix::det() const {
    SquareMatrix B = *this;
    double determinant = 1.0;
    int sign = 1;

    for (size_t step = 0; step < dimension - 1; ++step) {
        // Check if the pivot element B[step][step] is zero
        if (std::fabs(B(step, step)) < std::numeric_limits<double>::epsilon()) {
            // Try to find a non-zero element in the same column (row permutations)
            bool found = false;
            for (size_t row = step + 1; row < dimension; ++row) {
                if (std::fabs(B(row, step)) > std::numeric_limits<double>::epsilon()) {
                    B.permutateRows(step, row);
                    sign *= -1;
                    found = true;
                    break;
                }
            }
            if (!found) {
                // Try to find a non-zero element in the same row (column permutations)
                for (size_t col = step + 1; col < dimension; ++col) {
                    if (std::fabs(B(step, col)) > std::numeric_limits<double>::epsilon()) {
                        B.permutateCols(step, col);
                        sign *= -1;
                        found = true;
                        break;
                    }
                }
            }
            // If no non-zero element is found, the determinant is zero
            if (!found) {
                return 0.0;
            }
        }

        // Proceed with Gaussian elimination
        for (size_t row = step + 1; row < dimension; ++row) {
            double coeff = -B(row, step) / B(step, step);
            for (size_t col = step; col < dimension; ++col) {
                B(row, col) += B(step, col) * coeff;
            }
        }
    }

    // Calculate the product of the diagonal elements
    for (size_t i = 0; i < dimension; ++i) {
        determinant *= B(i, i);
    }

    return determinant * sign;
}*/

std::optional<size_t> SquareMatrix::find_zero() const
{
    std::optional<size_t> res{std::nullopt};

    for (size_t i = 0; i < dimension; ++i) {
        if (std::fabs(data[i * dimension + i]) < std::numeric_limits<double>::epsilon()) {
            res = i;
            break;
        }
    }
    return res;
}

SquareMatrix SquareMatrix::prepare_for_Gaus(double& sign) const
{
    SquareMatrix res = *this;
    struct Diag {
        size_t place;
        size_t direction;
        size_t idx_x;
        size_t idx_y;
    };

    std::vector<Diag> diag(dimension);
    size_t curElem = 0;
    sign = 1;

    while (!std::isnan(res[0][0])) {
        //std::cout << "Check diag:\n" << res << "\n";
        std::optional<size_t> zero_item = res.find_zero();
        if (!zero_item.has_value()) {
            break;
        }

        bool new_perm = true;
        Diag d = { zero_item.value(), 0, 0, 0 };
        do {
            while (d.direction == 0) {
                // check cols for permutation
                if (d.idx_x == d.place) {
                    ++d.idx_x;
                }
                if (d.idx_x == dimension) {
                    d.direction = 1;
                    d.idx_x = d.place;
                    d.idx_y = 0;
                    break;
                }
                // Probably exchange d.place column with d.idx_x column. Needs additional check
                if (std::fabs(res[d.idx_x][d.place]) < std::numeric_limits<double>::epsilon()) {
                    ++d.idx_x;
                    continue;
                }
                if ((d.idx_x < d.place) && std::fabs(res[d.place][d.idx_x]) < std::numeric_limits<double>::epsilon()) {
                    // If new diagonal value is zero. Needs to skip this case
                    ++d.idx_x;
                    continue;
                }
                break;
            }
            while (d.direction == 1) {
                // check rows for permutation
                if (d.idx_y == d.place) {
                    ++d.idx_y;
                }
                if (d.idx_y == dimension) {
                    // reject current permutation case
                    new_perm = false;
                    break;
                }
                // Probably exchange d.place row with d.idx_y row. Needs additional check
                if (std::fabs(res[d.place][d.idx_y]) < std::numeric_limits<double>::epsilon()) {
                    ++d.idx_y;
                    continue;
                }
                if ((d.idx_y < d.place) && std::fabs(res[d.idx_y][d.place]) < std::numeric_limits<double>::epsilon()) {
                    // If new diagonal value is zero. Needs to skip this case
                    ++d.idx_y;
                    continue;
                }
                break;
            }
            if (new_perm) {
                if (d.direction == 0) {
                    res.permutate_columns(d.place, d.idx_x);
                } else {
                    res.permutate_rows(d.place, d.idx_y);
                }
                //std::cout << "Permutated matrix:\n" << res << "\n";
                sign = sign * -1.;
                diag[curElem++] = d;
                break;
            } else {
                if (curElem) {
                    --curElem;
                    if (diag[curElem].direction == 0) {
                        res.permutate_columns(diag[curElem].place, diag[curElem].idx_x);
                    } else {
                        res.permutate_rows(diag[curElem].place, diag[curElem].idx_y);
                    }
                    //std::cout << "Reconstructed matrix:\n" << res << "\n";
                    sign = sign * -1.;
                    d = diag[curElem];
                } else {
                    // Correct permutation was not found
                    res[0][0] = NAN;
                    break;
                }
            }
        } while (true);
    }

    return res;
}

double SquareMatrix::det() const
{
    double d = 0.;
    do {
        // Check rows
        for (size_t i = 0; i < dimension; ++i) {
            bool zero_row = true;
            for (size_t j = 0; j < dimension; ++j) {
                if (std::fabs(data[i * dimension + j]) > std::numeric_limits<double>::epsilon()) {
                    zero_row = false;
                    break;
                }
            }
            if (zero_row) {
                break;
            }
        }

        // Check columns
        for (size_t i = 0; i < dimension; ++i) {
            bool zero_col = true;
            for (size_t j = 0; j < dimension; ++j) {
                if (std::fabs(data[i + j * dimension]) > std::numeric_limits<double>::epsilon()) {
                    zero_col = false;
                    break;
                }
            }
            if (zero_col) {
                break;
            }
        }

        // Check diagonal and permutate if it needs
        double sign = 1.;
        SquareMatrix B = prepare_for_Gaus(sign);

        if (!std::isnan(B[0][0])) {
            //std::cout << "Transponded matrix:\n" << B << "\n";
            for (size_t step = 0; step < dimension - 1; step++) {
                for (size_t row = step + 1; row < dimension - 1; row++) {
                    double coeff = -B[row][step] / B[step][step];
                    for (size_t col = step; col < dimension; col++) {
                        B[row][col] += B[step][col] * coeff;
                    }
                }
            }
            //std::cout << "Converted matrix:\n" << B << "\n";
            d = sign;
            for (size_t i = 0; i < dimension; i++) {
                d *= B[i][i];
            }
        }
    } while (false);
    //std::cout << "--------========== Determinant is: " << d << " ==========--------\n\n";
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

void SquareMatrix::permutate_rows(size_t r1, size_t r2)
{
    if (r1 >= dimension || r2 >= dimension) {
        assert(false);
        return;
    }
    for (size_t i = 0; i < dimension; ++i) {
        double d = data[i + r1 * dimension];
        data[i + r1 * dimension] = data[i + r2 * dimension];
        data[i + r2 * dimension] = d;
    }
}

void SquareMatrix::permutate_columns(size_t c1, size_t c2)
{
    if (c1 >= dimension || c2 >= dimension) {
        assert(false);
        return;
    }
    for (size_t i = 0; i < dimension; ++i) {
        double d = data[c1 + i * dimension];
        data[c1 + i * dimension] = data[c2 + i * dimension];
        data[c2 + i * dimension] = d;
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
    return !(std::fabs(det()) > std::numeric_limits<double>::epsilon());
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
