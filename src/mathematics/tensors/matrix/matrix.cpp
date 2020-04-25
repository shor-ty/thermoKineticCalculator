/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     T hermo       | Open Source Thermo-Kinetic Library
  c-o-o-c     K iknetic     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Constructor.

    TKC is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or 
    (at your option) any later version.

    TKC is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TKC; if not, see <http://www.gnu.org/licenses/>

\*---------------------------------------------------------------------------*/

#include "definitions.hpp"
#include "matrix.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::Matrix::Matrix()
:
    Tensor(0, 0)
{}


TKC::Matrix::Matrix
(
    const size_t rows,
    const size_t cols,
    const scalar value
)
:
    Tensor(rows, cols, value)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::Matrix::~Matrix()
{}


// * * * * * * * * * * * * * * Operator Functions  * * * * * * * * * * * * * //


TKC::Matrix TKC::Matrix::operator*(const Matrix& B) const
{
    //- Calculation of the inner product of two matrices 
    //  The result is a matrix again C_ij = A_ik * B_kj
    
    //- Check if A cols == B rows
    if (cols() != B.rows())
    {
        ErrorMsg
        (
            "    The cols of matrix A does not match the rows of matrix B",
            __FILE__,
            __LINE__
        );
    }

    //- We are writing A * B it is non - cummutative
    //  The Matrix A is the object
    //  The Matrix B is the operators argument
    const Matrix& A = *this;

    //- The matrix multiplication returns a matrix C with A.rows & B.cols
    //  Index notation from above [Holzmann]
    Matrix C(rows(), B.cols());

    for (size_t i = 0; i < rows(); ++i)
    {
        for (size_t j = 0; j < B.cols(); ++j)
        {
            scalar sum{0};
            
            for (size_t k = 0; k < cols(); ++k)
            {
                sum += A(i,k) * B(k,j);
            }
            C(i,j,sum);
        }
    }

    return C;
}


TKC::Vector TKC::Matrix::operator*(const Vector& b) const
{
    //- Inner Product of Matrix and Vector
    
    //- Check if A cols == b rows
    if (cols() != b.rows())
    {
        ErrorMsg
        (
            "    The cols of matrix A does not match the rows of vector b",
            __FILE__,
            __LINE__
        );
    }

    //- The matrix A is the object itself
    const Matrix& A = *this;

    //- The inner product of a matrix and vector results in a vector
    //  See [Holzmann] for further information
    Vector c(b.rows());


    for(size_t i=0; i<rows(); i++)
    {
        for (size_t j=0; j<b.rows(); j++)
        {
            c(i) += A(i,j) * b(j);
        }
    }

    return c;
}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * //

TKC::Matrix TKC::Matrix::T() const
{
    //- Cols and rows of actual matrix
    const size_t& row = this->rows();
    const size_t& col = this->cols();

    //- Transposed Matrix
    Matrix AT(col, row);

    //- Reference to actual matrix A
    const Matrix& A = *this;

    for (size_t r = 0; r < row; ++r)
    {
        for (size_t c = 0; c < col; ++c)
        {
            AT(c, r, A(r, c));
        }
    }

    return AT;

}


TKC::Matrix TKC::Matrix::inverse() const
{
    //- rows and cols of actual matrix
    const size_t row = rows();
    const size_t col = cols();

    //- Inverse only for squared matrix
    if (row != col)
    {
        ErrorMsg
        (
            "    You try to use A.inverse() for a non-squared matrix.\n"
            "    The inverse function is only available for squared matrices",
            __FILE__,
            __LINE__
        );
    }

    //- Copy of the actual matrix
    Matrix A = *this;

    //- GAUSS ELIMINATION

    //- First we need the identiy to build the inverse
    Matrix inverse = I();

    //- Make the Lower to zero
    for (unsigned int r = 0; r < row; r++)
    {
        //- Get first element
        const scalar element = A(r,r);        

        //- Divide the whole row by this element to get ONE at the element pos
        for (unsigned int c = 0; c < col; c++)
        {
            const scalar oldElem = A(r, c);
            const scalar oldElemI = inverse(r, c);

            A(r, c, oldElem / element); 
            inverse(r, c, oldElemI / element);
        }

        //- Now substract from the others to generate ZERO below the element
        for (unsigned int rr = r+1; rr < row; rr++)
        {
            const scalar multiplicator = A(rr, r) / A(r, r);

            for (unsigned int cc = 0; cc < col; cc++)
            {
                const scalar actualElement = A(rr, cc);
                const scalar modifiedElement = A(r, cc);

                const scalar actualElementI = inverse(rr, cc);
                const scalar modifiedElementI = inverse(r, cc);

                A(rr, cc, actualElement - modifiedElement * multiplicator);
                inverse
                (
                    rr, cc, actualElementI - modifiedElementI * multiplicator
                );
            }
        }
    }

    //- Make the Upper to zero. Diagonals are 1 now
    for (int r = row-1; r >= 0; r--)
    {
        //- Now substract from the others to generate ZERO above the element
        for (int rr = r-1; rr >= 0; rr--)
        {
            const scalar multiplicator = A(rr, r);

            for (int cc = col-1; cc >= 0; cc--)
            {
                const scalar actualElement = A(rr, cc);
                const scalar modifiedElement = A(r,cc);

                const scalar actualElementI = inverse(rr, cc);
                const scalar modifiedElementI = inverse(r,cc);

                A(rr, cc, actualElement - modifiedElement * multiplicator);
                inverse
                (
                    rr, cc, actualElementI - modifiedElementI * multiplicator
                );
            }
        }
    }

    return inverse;
}



TKC::Matrix TKC::Matrix::I() const
{
    size_t n{0};

    if (rows() > cols())
    {
        n = rows();
    }
    else
    {
        n = cols();
    }

    //- Build n x n matrix with zero
    Matrix tmp(n, n);

    for (unsigned int i = 0; i < n; i++)
    {
        tmp(i, i, 1);
    }

    return tmp;
}


TKC::Matrix TKC::Matrix::I(const size_t n) const
{
    //- Build n x n matrix with zero
    Matrix tmp(n, n);

    for (unsigned int i = 0; i < n; ++i)
    {
        tmp(i, i, 1);
    }

    return tmp;
}


TKC::Matrix TKC::Matrix::LT() const
{
    if (this->rows() != this->cols())
    {
        ErrorMsg
        (
            "Matrix is not a squared matrix",
            __FILE__,
            __LINE__
        );
    }

    Matrix L = this->I();

    size_t n = this->rows();

    //- Fill L with values of matrix
    for (size_t i=0; i<n; i++)
    {
        for (size_t j=0; j<i; j++)
        {
            L(i,j) = this->operator()(i,j);
        }
    }

    return L;
}


TKC::Matrix TKC::Matrix::UT() const
{
    if (this->rows() != this->cols())
    {
        ErrorMsg
        (
            "Matrix is not a squared matrix",
            __FILE__,
            __LINE__
        );
    }

    Matrix U = this->I();

    size_t n = this->rows();

    //- Fill U with values of matrix
    for (size_t i=0; i<n; i++)
    {
        for (size_t j=i; j<n; j++)
        {
            U(i,j) = this->operator()(i,j);
        }
    }

    return U;
}


void TKC::Matrix::LU() const
{
    //- Temporary matrix U0
    Matrix U0 = *this;

    //- Temporary matrix L0
    Matrix L0 = this->I();


    U0();

    L0();
}

// * * * * * * * * * * * * Special Matrix Functions  * * * * * * * * * * * * //

/*void TKC::Matrix::polynomCoefficients
(
    const scalarField& x,
    const scalarField& y,
    scalarField& coeffs
) const
{
    int i{0};
    int j{0};
    int k{0};
    int n = x.size();

    scalar phi;
    scalar ff;
    scalar b;

    scalarField s(n, scalar(0));
    scalarField c(n, scalar(0));
    
    s[n-1] = -x[0];

    for (i = 1; i < n; i++)
    {
        for (j = n-1-i; j < n-1; j++)
        {
            s[j] -= x[i] * s[j+1];
        }

        s[n-1] -= x[i];
    }

    for (j = 0; j < n; j++)
    {
        phi = n;

        for (k = n-1; k > 0; k--)
        {
            phi = k * s[k] + x[j] * phi;
        }

        ff = y[j] / phi;

        b = scalar(1);

        for (k = n-1; k >= 0; k--)
        {
            c[k] += b*ff;
            b = s[k] + x[j] * b;
        }
    }

    forAll(c, coeff)
    {
        Info<< "Coeff: " << coeff << "\n";
    }
}*/

// ************************************************************************* //
