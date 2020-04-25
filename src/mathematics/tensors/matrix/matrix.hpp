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

Class
    TKC::Matrix
    
Description
    Abstract TKC::Matrix class for building matrix and dealing with them
    The matrix of rank 0 are scalars
    The matrix of rank 1 are vectors
    The matrix of rank 2 are matrices (n x m)


SourceFiles
    matrix.cpp

\*---------------------------------------------------------------------------*/

#ifndef Matrix_hpp
#define Matrix_hpp

#include "definitions.hpp"
#include "tensor.hpp"
#include "vector.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class Vector;


/*---------------------------------------------------------------------------*\
                            Class Matrix Declaration
\*---------------------------------------------------------------------------*/

class Matrix
:
    public Tensor
{
    public:

        //- Constructor 
        Matrix();

        //- Constructor that creates a matrix n x m and init with the given
        //  value (default value = 0)
        Matrix(const size_t, const size_t, const scalar value = 0);

        //- Destructor
        ~Matrix();

        // Operator Functions


        // Arithmetric 

            //- Matrix multiplication (n x m) * (m x n) 
            //  The inner product A bullet B
            Matrix operator*(const Matrix&) const;

            //- Inner product of Matrix and Vector
            Vector operator*(const Vector&) const;


        // Member Functions

            //- Transpone the matrix A -> A^T
            Matrix T() const;

            //- Calcualte and return the inverse matrix of a squared matrix
            Matrix inverse() const;

            //- Return the identity matrix I with size of the call matrix
            Matrix I() const;

            //- Return the identity matrix I with given size (n x n)
            Matrix I(const size_t) const;

            //- Return Lower Triagonal of the matrix (LT has 1 at diag)
            //  Only for squared matrices
            Matrix LT() const;

            //- Return Upper Triagonal of the matrix
            //  Only for squared matrices
            Matrix UT() const;


            //- LU decomposition; decompose the Matrix A into L and U
            void LU() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Matrix_hpp included

// ************************************************************************* //
