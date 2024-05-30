#ifndef MATRIXINVERSE_H
#define MATRIXINVERSE_H

#include <vector>

struct Results;
class SquareMatrix;

class CBaseMatrixInverse
{
public:
    CBaseMatrixInverse();
    CBaseMatrixInverse(const CBaseMatrixInverse&);
    virtual ~CBaseMatrixInverse();
    virtual Results inverse(const SquareMatrix&, double precision = 0.) = 0;
};

class CLUPMatrixInverse : public CBaseMatrixInverse
{
public:
    CLUPMatrixInverse();
    CLUPMatrixInverse(const CLUPMatrixInverse&);
    virtual ~CLUPMatrixInverse();
    virtual Results inverse(const SquareMatrix&, double precision = 0.) override;
private:
    Results decompose(const SquareMatrix&);
    void LU_solve(Results& res);
};

class CShultsMatrixInverse : public CBaseMatrixInverse
{
public:
    CShultsMatrixInverse();
    CShultsMatrixInverse(const CShultsMatrixInverse&);
    virtual ~CShultsMatrixInverse();
    virtual Results inverse(const SquareMatrix&, double precision = 0.) override;
};

#endif // MATRIXINVERSE_H
