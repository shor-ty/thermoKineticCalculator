/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
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

Class
    AFC::LUDecompose
    
Description
    Decompose the matrix in a Lower (L) and Upper (U) matrix in order to solve
    it with forward and backward substitution. A = LU 
    Decomposition is done with pivoting. The matrix has to be a square matrix

    William H. Press, Saul A. Teukolsky, William T. Vetterling and Brian P.
    Flannery "Numerical Recipes - The Art of Scientif Computing, Third Edition"

SourceFiles
    LUDecompose.cpp

\*---------------------------------------------------------------------------*/

#ifndef LUDecompose_hpp
#define LUDecompose_hpp

#include "definitions.hpp"
#include "matrix.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class LUDecompose Declaration
\*---------------------------------------------------------------------------*/

class LUDecompose
{
    private:

        //- Rows of matrix 
        size_t i_;

        //- Cols of matrix
        size_t j_;

        //- Size of the matrix (Rows and cols has to be equal)
        size_t n_;

        //- Original matrix (needed for iterative improvment)
        const Matrix A_;

        //- LU matrix (store the decomposition)
        Matrix LU_;

        //- Stores the permutation for solution vector x and source b
        Vector permut_;

        //- Parity used for determinant
        scalar parity_{scalar(1)};



    public:

        //- Constructor 
        LUDecompose(Matrix&);

        //- Destructor
        ~LUDecompose();


        // Operator Functions

            //- Output the LU matrix
            void operator()() const;


        // Member functions

            //- Solve the linear set of equations A x = b
            void solve(Vector&, Vector&);

            //- Solution improve. The solution optained by solve() can be
            //  improved iteratively
            void improveSolution(const Vector&, const Vector&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // LUDecompose_hpp included

// ************************************************************************* //
