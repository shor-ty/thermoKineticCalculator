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
    TKC::Dimensioned
    
Description
    Abstract TKC::Dimensioned class for checking the dimensions of quantities


SourceFiles
    dimensioned.cpp

\*---------------------------------------------------------------------------*/

#ifndef Dimensioned_hpp
#define Dimensioned_hpp

#include "scalar.hpp"
#include "vector.hpp"
#include "matrix.hpp"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                            Class Dimensioned Declaration
\*---------------------------------------------------------------------------*/

template<typename Type>
class Dimensioned
{
    private:

        // Debug swithc
        bool debug_{true};

        //- The dimension field
        //  SI units
        //  1: [kg]
        //  2: [m]
        //  3: [s]
        //  4: [K]
        scalarField dimension_;

        //- Data name
        word name_;

        //- Data value
        Type value_;


    public:
        
        //- Constructor with dimensions
        Dimensioned(size_t, size_t, size_t, size_t);   

        //- Constructor with dimensions
        Dimensioned(const scalarField&);

        //- Destructor
        ~Dimensioned();

        // Operator Functions


        // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Dimensioned_hpp included

// ************************************************************************* //
