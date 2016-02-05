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

\*---------------------------------------------------------------------------*/

#include "chemistry.hpp"
#include "thermo.hpp"
#include "constants.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Chemistry::Chemistry
(
    const string& fileName 
)
{
    Info<< "Constructor Chemistry\n" << endl;

    ChemistryReader chemReader(fileName);
    
    chemReader.read(chemData_);

    //- Reactions that include species i
    createSpeciesInReaction();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Chemistry::~Chemistry()
{
    Info<< "Destructor Chemistry\n" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool AFC::Chemistry::thermo()
{
    return (chemData_.thermo());
}


// * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * * //

AFC::scalar AFC::Chemistry::calculateOmega
(
    const word& species,
    const scalar& T,
    map<word, scalar>& con,
    const Thermo& thermo
)
{
    //- Calculate source term omega
    return chemCalc_.calculateOmega(species, T, con, thermo, chemData_);
}


// * * * * * * * * * * * * * * * Create Functions  * * * * * * * * * * * * * //

void AFC::Chemistry::createSpeciesInReaction()
{
    //- Species list
    const wordList& species = chemData_.species();

    //- Reaction no.
    const int& nReac = chemData_.nReac();

    //- Wordmatrix that contains all species in each reaction
    const wordMatrix& speciesInReaction = chemData_.speciesInReactions();

    forAll(species, s)
    {
        //- If found species in reaction -> true
        bool found{false};

        //- Loop through all elementar reactions 
        for(int r=0; r<nReac; r++)
        {
            //- Loop through the species i in elementar reaction r
            for(unsigned int i=0; i<speciesInReaction[r].size(); i++)
            {
                if
                (
                  ! found
                 && speciesInReaction[r][i] == species[s]
                )
                {
                    found = true;

                    //- Insert reaction no to matrix
                    chemData_.setReacNumbers(species[s], r);
                }
            }

            //- Reset
            found = false;
        }
    }
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

AFC::wordList AFC::Chemistry::species() const
{
    return chemData_.species();
}
/*

int AFC::Chemistry::nReac() const
{
    return chemData_.nReac();
}


AFC::word AFC::Chemistry::elementarReaction
(
    const int& r
) const
{
    return chemData_.elementarReaction(r);
}


AFC::scalarField AFC::Chemistry::reacNoForSpecies
(
    const int& s
) const
{
    return chemData_.reacNoForSpecies(s);
}


AFC::scalarField AFC::Chemistry::k() const
{
    return chemData_.k();
}*/


/*AFC::wordMatrix AFC::Chemistry::speciesInReactions() const
{
    return chemData_.speciesInReactions();
}*/


/*AFC::scalar AFC::Chemistry::k
(
    const int& reacNo 
) const
{
    return chemData_.k(reacNo);
}*/


// ************************************************************************* //
