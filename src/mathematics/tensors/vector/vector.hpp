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
    TKC::Vector
    
Description
    Abstract TKC::Vector class that inherits all vector calculations 


SourceFiles
    vector.cpp

\*---------------------------------------------------------------------------*/

#ifndef Vector_hpp
#define Vector_hpp

#include "definitions.hpp"
#include "matrix.hpp"
#include "tensor.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class Matrix;


/*---------------------------------------------------------------------------*\
                            Class Vector Declaration
\*---------------------------------------------------------------------------*/

class Vector
:
    public Tensor
{
    public:

        //- Constructor 
        Vector();

        //- Copy constructor
        Vector(const Vector&);

        //- Constructor that creates a row vector n x 1 and init with the 
        //  given value.
        Vector(const size_t, const scalar value = 0);

        //- Constructor that creates a row vector n x 1 with an scalarField
        Vector(const scalarField&);

        //- Destructor
        ~Vector();


        //- Arithmetic

            //- Inner Product of Vector and Matrix
            Vector operator*(const Matrix&) const;

            //- Subtract one vector to the actual one
            void operator-=(const Vector&);

            //- Subtract one vector to the actual one and return the new vector
            const Vector operator-(const Vector&) const;


        // Operator Functions
        
            //- Return the value at position
            scalar operator()(const size_t) const;
        
            //- Set the value at position
            scalar& operator()(const size_t);

            //- Insert the value at position
            void operator()(const size_t, const scalar);

            //- Show the vector
            void operator()() const;

            //- Arithmetic += (const Matrix&)
            /*Vector operator()= 
            (
                const Matrix&
            );*/


        // Member Functions

            //- Return the size of the vector
            size_t size() const;

            //- Return the values of the vector
            scalarField values() const;


        // Calculation Functions

            //- Transpose the vector b -> b^T
            Vector T() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Vector_hpp included

// ************************************************************************* //
