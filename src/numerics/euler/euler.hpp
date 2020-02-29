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
    AFC::Euler
    
Description
    Abstract AFC::Euler class for building and calculating matrices

SourceFiles
    numerics.cpp

\*---------------------------------------------------------------------------*/

#ifndef Euler_hpp
#define Euler_hpp

#include "typedef.hpp"
#include "matrix.hpp"
//#include "jacobian.hpp"
#include "stepStatus.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Euler Declaration
\*---------------------------------------------------------------------------*/

class Euler
{
    private:

        // Debug switch
        bool debug_{false};

        //- Error of each species
        mutable map<word, scalar> err_;


    public:

        //- Constructor 
        Euler(const size_t);

        //- Destructor
        ~Euler();


        // Member functions

            //- Solve using euler algorithm
            scalar solve
            (
                const scalar,
                const scalar,
                const map<word, scalar>&,
                const map<word, scalar>&,
                map<word, scalar>&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Euler_hpp included

// ************************************************************************* //
