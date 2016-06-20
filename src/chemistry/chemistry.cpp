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
    ChemistryReader chemReader(fileName);
    
    chemReader.read(chemData_);

    //- Reactions that include species i
    createSpeciesInReaction();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Chemistry::~Chemistry()
{
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


void AFC::Chemistry::calculateKf
(
    const int& r,
    const scalar& T
)
{
    chemCalc_.calculateKf(r, T, chemData_);
}


void AFC::Chemistry::calculateKc
(
    const int& r,
    const scalar& T,
    const Thermo& thermo
)
{
    chemCalc_.calculateKc(r, T, thermo, chemData_);
}


void AFC::Chemistry::calculateKb()
{
    chemCalc_.calculateKb(chemData_);
}


// * * * * * * * * * * * * * * * Update Functions  * * * * * * * * * * * * * //

void AFC::Chemistry::updateM
(
    const scalar& M
)
{
    chemData_.updateM(M);
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

AFC::wordList AFC::Chemistry::elements() const
{
    return chemData_.elements();
}


AFC::wordList AFC::Chemistry::species() const
{
    return chemData_.species();
}


AFC::scalar AFC::Chemistry::kf() const
{
    return chemData_.kf();
}


AFC::scalar AFC::Chemistry::kb() const
{
    return chemData_.kb();
}


AFC::scalar AFC::Chemistry::Kc() const
{
    return chemData_.Kc();
}


AFC::word AFC::Chemistry::elementarReaction
(
    const int& r
) const
{
    return chemData_.elementarReaction(r);
}


AFC::stringList AFC::Chemistry::elementarReaction() const
{
    return chemData_.elementarReaction();
}


AFC::scalarList AFC::Chemistry::arrheniusCoeffs
(
    const int& reacNo
) const
{
    return chemData_.arrheniusCoeffs(reacNo);
}


AFC::scalarList AFC::Chemistry::LOWCoeffs
(
    const int& reacNo
) const
{
    return chemData_.LOWCoeffs(reacNo);
}


AFC::scalarList AFC::Chemistry::TROECoeffs
(
    const int& reacNo
) const
{
    return chemData_.TROECoeffs(reacNo);
}


AFC::scalarList AFC::Chemistry::SRICoeffs
(
    const int& reacNo
) const
{
    return chemData_.SRICoeffs(reacNo);
}


bool AFC::Chemistry::TBR
(
    const int& reacNo
) const
{
    return chemData_.TBR(reacNo);
}


bool AFC::Chemistry::ENHANCED
(
    const int& reacNo
) const
{
    return chemData_.ENHANCED(reacNo);
}


bool AFC::Chemistry::LOW
(
    const int& reacNo
) const
{
    return chemData_.LOW(reacNo);
}


bool AFC::Chemistry::TROE
(
    const int& reacNo
) const
{
    return chemData_.TROE(reacNo);
}


bool AFC::Chemistry::SRI
(
    const int& reacNo
) const
{
    return chemData_.SRI(reacNo);
}


AFC::map<AFC::word, AFC::scalar> AFC::Chemistry::enhancedFactors
(
    const int& reacNo
) const
{
    return chemData_.enhancedFactors(reacNo);
}


AFC::wordList AFC::Chemistry::enhancedSpecies
(
    const int& reacNo
) const
{
    return chemData_.enhancedSpecies(reacNo);
}


AFC::scalar AFC::Chemistry::dH() const
{
    return chemData_.dH();
}


AFC::scalar AFC::Chemistry::dS() const
{
    return chemData_.dS();
}


AFC::scalar AFC::Chemistry::dG() const
{
    return chemData_.dG();
}


// ************************************************************************* //
