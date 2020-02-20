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

#include "math.h" 
#include "typedef.hpp"
#include "LUDecompose.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::LUDecompose::LUDecompose(Matrix& A)
:
    i_(A.rows()),

    j_(A.cols()),

    A_(A),

    LU_(A),

    permut_(i_, 0)
{
    //- Check for square matrix
    if (i_ != j_)
    {
        FatalError 
        (
            "*** Error in LUDecompose::LUDecompose - Not a square matrix",
            __FILE__,
            __LINE__
        );
    }
    else
    {
        n_ = i_;
    }

    //- Temporary variables
    size_t i{0}, j{0}, k{0}, iMax{0};

    //- Temporary variable to check the pivot element and for switching
    scalar temp{0}, large{0}, pivot{0};

    //- Vector field that stores the implicit scaling of each row
    Vector scale(n_, 0);

    //- Search for the largest value in each row for scaling the matrix
    //  i := col
    //  j := row
    for (i=0; i<n_; i++)
    {
        //- Set to zero for new search
        large = scalar(0); 

        for (j=0; j<n_; j++)
        {
            //- Checking always absolute values
            if (fabs(A(i,j)) > large)
            {
                large = LU_(i,j);
            }
        }

        //- If in one row we have only zeros
        if (large == scalar(0))
        {
            FatalError 
            (
                "*** Error in LUDecompose::LUDecompose - Singular matrix",
                __FILE__,
                __LINE__
            );
        }

        //- Set scale factor for ith col 
        scale(i) = scalar(1)/large;
    }

    //- Checkout where the pivot element is
    for (k=0; k<n_; k++)
    {
        //- Reset
        large = scalar(0); 

        iMax = k;

        for (i=k; i<n_; i++)
        {
            temp = scale(i) * fabs(LU_(i,k));

            if (temp > large)
            {
                large = temp;

                //- Save position (row) of pivot 
                iMax = i;
            }  
        }

        //- Shift the entries if the pivot element is not already in the diag
        if (k != iMax)
        {
            for (j=0; j<n_; j++)
            {
                temp = LU_(iMax,j);

                //- Put all entries of kth row into actual pivot row
                LU_(iMax,j) = LU_(k,j);

                //- Put all entries of pivot row into actual kth row
                LU_(k,j) = temp;
            }

            //- Change pairity 
            parity_ *= scalar(-1);

            //- Shift also the scale factors
            temp = scale(iMax);

            scale(iMax) = scale(k);

            scale(k) = temp;
        }

        //- Store permutation
        permut_(k) = iMax;

        //- Just to make sure that the pivot element is not zero
        if (LU_(k,k) == scalar(0))
        {
            //- TODO 
            LU_(k,k) = 1e-40; //SMALL;
        }

        //- Build the LU matrix
        for (i=k+1; i<n_; i++)
        {
            temp = LU_(i,k) /= LU_(k,k);

            //- Reduce remaining submatrix
            for (j=k+1; j<n_; j++)
            {
                LU_(i,j) -= temp * LU_(k,j);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::LUDecompose::~LUDecompose()
{}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

void AFC::LUDecompose::solve(Vector& b, Vector& x)
{
    //- Temporary variables
    size_t i{0}, ii{0}, ip{0}, j{0};

    scalar sum{0};

    //- First check size of vectors compared to matrix
    if (x.size() != n_)
    {
        FatalError
        (
            "*** Error in LUDecompose::solve - Matrix size does not match"
            " with the vector x",
            __FILE__,
            __LINE__
        );
    }
    else if (b.size() != n_)
    {
        FatalError
        (
            "Matrix size does not match with vector b",
            __FILE__,
            __LINE__
        );
    }

    //- Setting x values equal to b 
    for (i=0; i<n_; i++)
    {
        x(i) = b(i);
    }

    //- Forward substitution
    for (i=0; i<n_; i++)
    {
        ip = permut_(i);

        sum = x(ip);

        x(ip) = x(i);

        if (ii != 0)
        {
            for (j=ii-1; j<i; j++)
            {
                sum -= LU_(i,j) * x(j);
            }
        }
        else if (sum != scalar(0))
        {
           ii = i+1; 
        }

        x(i) = sum;
    }

    //- Backward substitution 
    for (int i=n_-1; i>=0; i--)
    {
        sum = x(i);

        for (j=i+1; j<n_; j++)
        {
            sum -= LU_(i,j) * x(j);
        }

        //- Store the component of the solution vector
        x(i) = sum / LU_(i,i);
    }
}


void AFC::LUDecompose::improveSolution(Vector& b, Vector& x)
{
    //- Solve the system Ax to get b* (x is the calculated solution)
    Vector bstar = A_ * x;

    //- The error is simply the difference between the real solution and bstar
    bstar -= b;

    //- Solve for the error term
    solve(bstar, bstar);
    
    //- Substract the error from the solution
    x-= bstar;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
