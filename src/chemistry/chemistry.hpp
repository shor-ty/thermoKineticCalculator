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
    AFC::Chemistry
    
Description
    Abstract AFC::Chemistry class for chemistry data and calculation

SourceFiles
    chemistry.cpp

\*---------------------------------------------------------------------------*/

#ifndef Chemistry_hpp
#define Chemistry_hpp

#include "chemistryReader.hpp"
#include "chemistryData.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Chemistry Declaration
\*---------------------------------------------------------------------------*/

class Chemistry
{
    private:

        // Private pointer data

            //- chemistryData obj
            ChemistryData chemData_;
    

    public:

        //- Constructor
        Chemistry
        (
            const string&
        );

        //- Destructor
        ~Chemistry();


        // Member Functions
            
            //- Read chemistry file
            void readChemistry();

            //- Return chemistryData::thermo_
            bool thermo();


        // Return functions
        
            //- Retrun species
            wordList species() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Chemistry_hpp included

// ************************************************************************* //
