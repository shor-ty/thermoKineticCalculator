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
    TKC::Scalar
    
Description
    Abstract TKC::Scalar class that inherits all scalar calculations 


SourceFiles
    scalar.cpp

\*---------------------------------------------------------------------------*/

#ifndef Scalar_hpp
#define Scalar_hpp

#include "definitions.hpp"
#include "matrix.hpp"
#include "tensor.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class Matrix;


/*---------------------------------------------------------------------------*\
                            Class Scalar Declaration
\*---------------------------------------------------------------------------*/

class Scalar
:
    public Tensor
{
    public:

        //- Constructor 
        Scalar();

        //- Copy constructor
        Scalar(const Scalar&);

        //- Constructor that creates a row scalar n x 1 and init with the 
        //  given value.
        Scalar(const size_t, const scalar value = 0);

        //- Constructor that creates a row scalar n x 1 with an scalarField
        Scalar(const scalarField&);

        //- Destructor
        ~Scalar();

        // Operator Functions

            Scalar operator*(const Matrix&) const;
        
            //- Return the value at position
            scalar operator()(const size_t) const;

            //- Insert the value at position
            void operator()(const size_t, const scalar);

            //- Show the scalar
            void operator()() const;

            //- Arithmetic += (const Matrix&)
            /*Scalar operator()= 
            (
                const Matrix&
            );*/


        // Member Functions

            //- Return the size of the scalar
            size_t size() const;

            //- Return the values of the scalar
            scalarField values() const;


        // Calculation Functions

            //- Transpose the scalar b -> b^T
            Scalar T() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Scalar_hpp included

// ************************************************************************* //
