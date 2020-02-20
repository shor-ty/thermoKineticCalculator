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
    AFC::Dimensioned
    
Description
    Abstract AFC::Dimensioned class for checking the dimensions of quantities


SourceFiles
    dimensioned.cpp

\*---------------------------------------------------------------------------*/

#ifndef Dimensioned_hpp
#define Dimensioned_hpp

#include "scalar.hpp"
#include "vector.hpp"
#include "matrix.hpp"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
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

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Dimensioned_hpp included

// ************************************************************************* //
