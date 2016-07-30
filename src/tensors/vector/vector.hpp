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

Class
    AFC::Vector
    
Description
    Abstract AFC::Vector class that inherits all vector calculations 


SourceFiles
    vector.cpp

\*---------------------------------------------------------------------------*/

#ifndef Vector_hpp
#define Vector_hpp

#include "typedef.hpp"
#include "matrix.hpp"
#include "tensor.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
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
    private:

        // Debug swithc
        bool debug_{true};

        // Private Data


    public:

        //- Constructor 
        Vector();

        //- Copy constructor
        Vector
        (
            const Vector&
        );

        //- Constructor that creates a row vector n x 1 and init with the 
        //  given value.
        Vector
        (
            const size_t,
            const scalar value = 0
        );

        //- Constructor that creates a row vector n x 1 with an scalarField
        Vector
        (
            const scalarField&
        );

        //- Destructor
        ~Vector();

        // Operator Functions

            Vector operator*
            (
                const Matrix&
            ) const;
        
            //- Return the value at position
            scalar operator()
            (
                const size_t&
            ) const;

            //- Insert the value at position
            void operator()
            (
                const size_t,
                const scalar
            );

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


        // Calculation Functions

            //- Transpose the vector b -> b^T
            Vector T() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Vector_hpp included

// ************************************************************************* //
