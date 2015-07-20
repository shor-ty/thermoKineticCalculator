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
    AFC::MixtureFractionReader    

Description
    Reading the AFCDict file

SourceFiles
    mixtureFractionReader.cpp

\*---------------------------------------------------------------------------*/

#ifndef MixtureFractionReader_hpp
#define MixtureFractionReader_hpp

#include "stringManipulator.hpp"
#include "mixtureFraction.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{


class MixtureFraction;

/*---------------------------------------------------------------------------*\
                      Class MixtureFractionReader Declaration
\*---------------------------------------------------------------------------*/

class MixtureFractionReader
:
    public StringManipulator
{
    private:

        // Private data

            //- MixtureFraction file
            string file_;


    public:

        // Constructor and Destructor

            //- Constructor with file string and MixtureFraction:: obj adress
            MixtureFractionReader
            (
                const string&
            );

            //- Destructor
            ~MixtureFractionReader();

        
        // Member functions

            //- Read mixtureFraction file and delegate data
            void read
            (
                MixtureFraction&
            );


        // Helper functions
        
            //- Find line number of keyword

            //- Return string between '/' and '/'


        // Data manipulation functions
        
            //- Manipulate elementar reaction string and analyze reaction
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
