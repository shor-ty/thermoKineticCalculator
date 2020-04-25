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
    TKC::Tensor
    
Description
    Abstract TKC::Tensor class for building tensor and dealing with them
    The tensor of rank 0 are scalars
    The tensor of rank 1 are vectors
    The tensor of rank 2 are matrices (n x m)


SourceFiles
    tensor.cpp

\*---------------------------------------------------------------------------*/

#ifndef Tensor_hpp
#define Tensor_hpp

#include "definitions.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                            Class Tensor Declaration
\*---------------------------------------------------------------------------*/

class Tensor
{
    private:

        // Private Data

            //- Rows
            const size_t nRows_;

            //- Colums
            const size_t nCols_;

            //- Value
            scalarField values_;


    public:

        //- Constructor 
        Tensor();

        //- Constructor that create tensor n x m and init with zero
        Tensor(const size_t, const size_t, const scalar value = 0);

        //- Constructor that creates a tensor n x m using a scalarField
        Tensor(const size_t, const size_t, const scalarField&);

        //- Destructor
        ~Tensor();


        // Operator Functions

            //- Return the tensor element at the given positions (matrix)
            scalar operator()(size_t, size_t) const;

            //- Return the tensor element at the given positions (vector)
            scalar operator()(size_t) const;

            //- Assign operator
            scalar& operator()(const size_t, const size_t);

            //- Assign operator (b(i) = value)
            scalar& operator()(const size_t);

            //- Set tensor element at given position
            void operator()(size_t, size_t, scalar);

            //- Output the full tensor
            void operator()() const;


        // Member Functions

            //- Return cols j of the tensor (Spalten)
            size_t cols() const;

            //- Return rows i of the tensor (Zeilen)
            size_t rows() const;

            //- Return the values of the tensor
            scalarField values() const;

            //- Reset the tensor and set all values to zero
            void reset();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Tensor_hpp included

// ************************************************************************* //
