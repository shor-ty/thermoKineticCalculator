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
    AFC::Jacobian
    
Description
    Abstract AFC::Jacobian class for calculating the Jacobian matrix

SourceFiles
    jacobian.cpp

\*---------------------------------------------------------------------------*/

#ifndef Jacobian_hpp
#define Jacobian_hpp

#include "typedef.hpp"
#include "matrix.hpp"
#include "chemistry.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Jacobian Declaration
\*---------------------------------------------------------------------------*/

class Jacobian
{
    private:

        // Reference to the chemistry object
        const Chemistry& chem_;


    public:

        //- Constructor 
        Jacobian(const Chemistry&);

        //- Destructor
        ~Jacobian();


        // Member functions

            //- Calculate the Jacobian matrix
            void jacobian 
            (
                const scalar,
                const scalar,
                const scalar,
                const map<word, scalar>&,
                map<word, scalar>&,
                Matrix&
            );

            //- Derive the elementar reaction based on the species s and return
            //  the value of the derivation
            scalar derivationOfReaction
            (
                const word,
                const word,
                const wordList&,
                const wordList&,
                const map<word, int>&,
                const map<word, int>&,
                const scalar,
                const scalar,
                const map<word, scalar>&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Jacobian_hpp included

// ************************************************************************* //
