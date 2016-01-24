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
#include "chemistryCalc.hpp"
#include "thermo.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Chemistry Declaration
\*---------------------------------------------------------------------------*/

class Chemistry
{
    public:

        // Private reference data

            //- ChemistryData obj
            ChemistryData chemData_;

            //- ChemistryCalc obj
            ChemistryCalc chemCalc_;


        // Debug 
        bool debug{false};

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


        // Calculation Functions

            //- Update the reaction rate kf and kb
            void updatekfkb
            (
                const scalar&,
                const Thermo&
            );

            //- Calculate the source term of each species (omega)
            void omega
            (
                const map<word, scalar>&
            );


        // Create Functions

            //- Create field that contains all reaction no. for each species
            void createSpeciesInReaction();

            //- Insert reaction no. for species ::DELEGATED to CHEMISTRYDATA
            void insertReacNoForSpecies
            (
                const int&,
                const int&
            );


        // Return Functions
        
            //- Retrun species
            wordList species() const;

            //- Return no. of elementar reaction
            /*/int nReac() const;

            //- Return elementar reaction (as string) DELEGATED
            word elementarReaction
            (
                const int&
            ) const;

            //- Return all reaction no of species DELEGATED
            scalarField reacNoForSpecies
            (
                const int&
            ) const;

            //- Return reaction rates k
            scalarField k() const;

            //- Return wordMatrix for species in reactions
            wordMatrix speciesInReactions() const;

            //- Return reaction rate k of reaction no.
            scalar k
            (
                const int& 
            ) const;*/

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Chemistry_hpp included

// ************************************************************************* //
