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
#include "constants.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Chemistry::Chemistry
(
    const string& fileName 
)
{
    ChemistryReader chemReader(fileName);
    
    chemReader.read(chemData_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Chemistry::~Chemistry()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool AFC::Chemistry::thermo()
{
    return (chemData_.thermo());
}


// * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * * //

void AFC::Chemistry::k
(
    const scalar& T
)
{
    //- No. of reactions
    const int& reac = chemData_.nReac();

    scalarField k;

    //- Calculate reaction rates k
    for (int r=0; r<reac; r++)
    {
        //- Arrhenius coeffs
        const scalarField& arrCoeffs = chemData_.arrheniusCoeffs(r);

        //- Pre-exponential factor
        const scalar A = arrCoeffs[0];

        //- Expontent
        const scalar beta = arrCoeffs[1];

        //- Activation energy
        const scalar Ea  = arrCoeffs[2];  

        //- Normal calculation
        if
        (
            !chemData_.LOW(r)
         && !chemData_.TROE(r)
         && !chemData_.SRI(r)
         && !chemData_.ENHANCED(r)
         && !chemData_.TBR(r)
        )
        {
            //- ARRHENIUS EQUATION
            // * * * * * * * * * * * * * * * * * * //
            k.push_back(A * pow(T, beta) * exp(Ea/(AFC::Constants::R*T)));
            // * * * * * * * * * * * * * * * * * * //
            //
        }
    }

    //- Move calculated reaction rates into chemData_::reacRates_
    chemData_.update_k(k);
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

AFC::wordList AFC::Chemistry::species() const
{
    return chemData_.species();
}


int AFC::Chemistry::nReac() const
{
    return chemData_.nReac();
}


AFC::scalarField AFC::Chemistry::k() const
{
    return chemData_.k();
}


/*AFC::scalar AFC::Chemistry::k
(
    const int& reacNo 
) const
{
    return chemData_.k(reacNo);
}*/


// ************************************************************************* //
