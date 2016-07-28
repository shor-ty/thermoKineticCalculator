/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2015 Holzmann-cfd
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Constructor.

    AFC is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or 
    (at your option) any later version.

    AFC is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AFC; if not, see <http://www.gnu.org/licenses/>

\*---------------------------------------------------------------------------*/

#include "typedef.hpp"
#include "matrix.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Matrix::Matrix()
:
    Tensor(0, 0)
{
    if (debug_)
    {
        Info<< "Constructor Matrix()\n" << endl;
    }
}


AFC::Matrix::Matrix
(
    const size_t rows,
    const scalar value
)
:
    Tensor(rows, rows, value)
{
    if (debug_)
    {
        Info<< "Constructor Matrix (rows, rows, value)\n" << endl;
    }
}


AFC::Matrix::Matrix
(
    const size_t rows,
    const size_t cols,
    const scalar value
)
:
    Tensor(rows, cols, value)
{
    if (debug_)
    {
        Info<< "Constructor Matrix (rows, cols, value)\n" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Matrix::~Matrix()
{
    if (debug_)
    {
        Info<< "Destructor Matrix \n" << endl;
    }
}


// * * * * * * * * * * * * * * Operator Functions  * * * * * * * * * * * * * //

/*
AFC::scalarField AFC::Matrix::operator*
(
    const Matrix& B
) const
{
    //- Check if A cols == B rows
    if (nCols_ != B.rows())
    {
        FatalError
        (
            "    The cols of matrix A does not match the rows of matrix B",
            __FILE__,
            __LINE__
        );
    }

    //- The matrix before the * sign (A * B)
    const Matrix& A = *this;

    //- The matrix multiplication returns a vector b (scalarField)
    scalarField b(nRows_, 0);

    for (size_t r = 0; r < nRows_; r++)
    {
        for (size_t c = 0; c < nCols_; c++)
        {
            b[c] += A(r,c) * B(c,r);
        }
    }

    return b;
}*/


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * //

AFC::Matrix AFC::Matrix::T() const
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
/*AFC::Matrix AFC::Matrix::inverse() const
{
    //- Inverse only for squared matrix
    if (nRows_ != nCols_)
    {
        FatalError
        (
            "    The inverse function is only available for squared matrices",
            "\n    You try to use A.inverse() for a non-squared matrix. For\n"
            "    non-squared matrices you can use inverse"
            __FILE__,
            __LINE__
        );
    }

    //- Copy of the actual matrix
    Matrix A = *this;

    //- First we need the identiy to build the inverse
    Matrix inverse = I();

    //- Make the Lower to zero
    for (unsigned int r = 0; r < nRows_; r++)
    {
        //- Get first element
        const scalar element = A(r,r);        

        //- Divide the whole row by this element to get ONE at the element pos
        for (unsigned int c = 0; c < nCols_; c++)
        {
            const scalar oldElem = A(r, c);
            const scalar oldElemI = inverse(r, c);

            A(r, c, oldElem / element); 
            inverse(r, c, oldElemI / element);
        }

        //- Now substract from the others to generate ZERO below the element
        for (unsigned int rr = r+1; rr < nRows_; rr++)
        {
            const scalar multiplicator = A(rr, r) / A(r, r);

            for (unsigned int cc = 0; cc < nCols_; cc++)
            {
                const scalar actualElement = A(rr, cc);
                const scalar modifiedElement = A(r, cc);

                const scalar actualElementI = inverse(rr, cc);
                const scalar modifiedElementI = inverse(r, cc);

                A(rr, cc, actualElement - modifiedElement * multiplicator);
                inverse(rr, cc, actualElementI - modifiedElementI * multiplicator);
            }
        }
    }

    //- Make the Upper to zero. Diagonals are 1 now
    for (int r = nRows_-1; r >= 0; r--)
    {
        //- Now substract from the others to generate ZERO above the element
        for (int rr = r-1; rr >= 0; rr--)
        {
            const scalar multiplicator = A(rr, r);

            for (int cc = nCols_-1; cc >= 0; cc--)
            {
                const scalar actualElement = A(rr, cc);
                const scalar modifiedElement = A(r,cc);

                const scalar actualElementI = inverse(rr, cc);
                const scalar modifiedElementI = inverse(r,cc);

                A(rr, cc, actualElement - modifiedElement * multiplicator);
                inverse(rr, cc, actualElementI - modifiedElementI * multiplicator);
            }
        }
    }

    return inverse;
}


AFC::Matrix AFC::Matrix::I() const
{
    size_t n{0};

    if (nRows_ > nCols_)
    {
        n = nRows_;
    }
    else
    {
        n = nCols_;
    }

    //- Build n x n matrix with zero
    Matrix tmp(n, n);

    for (unsigned int i = 0; i < n; i++)
    {
        tmp(i, i, 1);
    }

    return tmp;
}


AFC::Matrix AFC::Matrix::I
(
    const size_t n 
) const
{
    //- Build n x n matrix with zero
    Matrix tmp(n, n);

    for (unsigned int i = 0; i < n; ++i)
    {
        tmp(i, i, 1);
    }

    return tmp;
}*/


// * * * * * * * * * * * * Special Matrix Functions  * * * * * * * * * * * * //

/*void AFC::Matrix::polynomCoefficients
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
