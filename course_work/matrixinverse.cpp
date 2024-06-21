#include <cmath>
#include "matrixinverse.h"
#include "matrix.h"
#include <QString>

QString matrix2str(const QString& caption, const SquareMatrix& m)
{
    QString result = caption + "\n";
    for (size_t j = 0; j < m.size(); ++j) {
        for (size_t i = 0; i < m.size(); ++i) {
            double d = m(i, j);
            if ((d > -1e5 && d < -1e-5) || (d > 1e-5 && d < 1e5) || d == 0) {
                result += QString().asprintf("%15.4lf", m[i][j]);
            } else {
                if ((d > -1e11 && d < -1e-11) || (d > 1e-11 && d < 1e11)) {
                    result += QString().asprintf("%15.4e", m[i][j]);
                } else {
                    result += QString().asprintf("%15.4lf", 0.);
                }
            }
        }
        result += "\n";
    }
    result += "\n";
    return result;
}

CBaseMatrixInverse::CBaseMatrixInverse()
{
}

CBaseMatrixInverse::CBaseMatrixInverse(const CBaseMatrixInverse&)
{
}

CBaseMatrixInverse::~CBaseMatrixInverse()
{
}

CLUPMatrixInverse::CLUPMatrixInverse() : CBaseMatrixInverse()
{
}

CLUPMatrixInverse::CLUPMatrixInverse(const CLUPMatrixInverse&) : CBaseMatrixInverse()
{
}

CLUPMatrixInverse::~CLUPMatrixInverse()
{
}

Results CLUPMatrixInverse::decompose(const SquareMatrix& matrix)
{
    Results res;
    res.complexity = 0;
    res.method = 0;
    res.determinant = matrix.det();
    res.LU = matrix;
    res.P = SquareMatrix::IdentityMatrix(matrix.size());
    for (int i = 0; i < matrix.size(); ++i) {
        double pivotValue = 0.;
        int pivot = -1;
        for (int row = i; row < matrix.size(); ++row) {
            if (fabs(res.LU[i][row]) > pivotValue) {
                pivotValue = fabs(res.LU[i][row]);
                pivot = row;
            }
            ++res.complexity;
        }
        if (pivotValue != 0) {
            res.P.permutate_rows(pivot, i);
            res.LU.permutate_rows(pivot, i);
            for (int j = i + 1; j < matrix.size(); ++j) {
                res.LU[i][j] /= res.LU[i][i];
                ++res.complexity;
                for (int k = i + 1; k < matrix.size(); ++k) {
                    res.LU[k][j] -= res.LU[i][j] * res.LU[k][i];
                    ++res.complexity;
                }
            }
        }
    }
    return res;
}

void CLUPMatrixInverse::LU_solve(Results& res)
{
    int n = static_cast<int>(res.LU.size());
    SquareMatrix b = *(const_cast<SquareMatrix*>(&(res.P))) * SquareMatrix::IdentityMatrix(n);
    for (int q = 0; q < n; q++) {
        for (int i = 0; i < n; i++) {           //forward substitution
            for (int j = 0; j < i; j++) {
                b[q][i] -= res.LU[j][i] * b[q][j];
                ++res.complexity;
            }
        }
        for (int i = n - 1; i > -1; i--) {      //back substitution
            res.X[i][q] = b[q][i];
            ++res.complexity;
            for (int j = i + 1; j < n; j++) {
                res.X[i][q] -= res.LU[j][i] * res.X[j][q];
                ++res.complexity;
            }
            res.X[i][q] = res.X[i][q] / res.LU[i][i];
            ++res.complexity;
        }
    }
    res.X.transpond();
}

Results CLUPMatrixInverse::inverse(const SquareMatrix& matrix, double)
{
    size_t n = matrix.size();

    Results res = decompose(matrix);

    res.X = res.LU;
    LU_solve(res);
    return res;
}

CShultsMatrixInverse::CShultsMatrixInverse() : CBaseMatrixInverse()
{
}

CShultsMatrixInverse::CShultsMatrixInverse(const CShultsMatrixInverse&) : CBaseMatrixInverse()
{

}

CShultsMatrixInverse::~CShultsMatrixInverse()
{
}

Results CShultsMatrixInverse::inverse(const SquareMatrix& matrix, double precision)
{
    Results res;
    res.complexity = 0;
    res.method = 1;
    res.precision = precision;
    res.determinant = matrix.det();
    SquareMatrix A0 = matrix;
    A0.transpond();
    A0.normalize();
    SquareMatrix E2 = SquareMatrix::IdentityMatrix(matrix.size());
    E2 = E2 * 2.;

    SquareMatrix inv = A0;
    int counter = 0;
    while (fabs( (matrix * inv).det() - 1.) >= precision) {       //пока |det(A * A[k](^-1)) - 1| >= EPS
        res.complexity += matrix.size() * matrix.size()* matrix.size();
        SquareMatrix prev = inv;            //A[k-1]
        inv = matrix * prev;                //A.(A[k-1]^(-1))
        res.complexity += matrix.size() * matrix.size()* matrix.size();
        inv = inv * (-1.);                  //-A.(A[k-1]^(-1))
        res.complexity += matrix.size() * matrix.size();
        inv = inv + E2;                     //2E - A.(A[k-1]^(-1))
        res.complexity += matrix.size() * matrix.size();
        inv = prev * inv;                   //(A[k-1]^(-1)).(2E - A.(A[k-1]^(-1)))
        res.complexity += matrix.size() * matrix.size()* matrix.size();
        ++counter;
        if (counter > 10000) {
            inv[0][0] = NAN;
            break;
        }
    }
    res.complexity += matrix.size() * matrix.size()* matrix.size();
    res.X = inv;
    res.LU[0][0] = counter;
    return res;
}
