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
    AFC::Dimensions
    
Description
    Abstract AFC::Dimensions class for checking the dimensions of quantities


SourceFiles
    dimensions.cpp

\*---------------------------------------------------------------------------*/

#ifndef Dimensions_hpp
#define Dimensions_hpp

#include "typedef.hpp"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                            Class Dimensions Declaration
\*---------------------------------------------------------------------------*/

class Dimensions
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
        const scalarField dimension_;


    public:

        //- Constructor with default dimensions [-] no unit
        Dimensions();

        //- Constructor with dimensions
        Dimensions
        (
            const scalarField&  
        );

        //- Destructor
        ~Dimensions();

        // Operator Functions


        // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Dimensions_hpp included

// ************************************************************************* //
