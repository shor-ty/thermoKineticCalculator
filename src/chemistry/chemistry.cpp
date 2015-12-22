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
    const scalar& T,
    const map<word, scalar>& speciesMol
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

        //- Pre-exponential factor [cm^3/mol/s]
        //  combustion.berkeley.edu
        const scalar A = arrCoeffs[0];

        //- Expontent
        const scalar beta = arrCoeffs[1];

        //- Activation energy [cal/mol]
        const scalar Ea  = arrCoeffs[2];  

        //- Calculate [M] if necessary
        
        scalar M{0};

        //- Calculate [M] out of literature [10] 
        //  TODO | check this
        if (!chemData_.ENHANCED(r))
        {
            M = calcM_Warnatz(speciesMol);    
        }
        else
        {
            M = calcM(speciesMol, r);    
        }
        

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
            k.push_back(arrhenius(A, beta, Ea, T));
        }
        //- LOW calculation
        /*else if
        (
            chemData_.LOW(r)
         && !chemData_.TROE(r)
         && !chemData_.SRI(r)
         && !chemData_.ENHANCED(r)
         && !chemData_.TBR(r)
        )
        {
            
        }
        //- TROE calculation
        else if
        (
            chemData_.LOW(r)
         && chemData_.TROE(r)
         && !chemData_.SRI(r)
         && !chemData_.ENHANCED(r)
         && !chemData_.TBR(r)
        )
        {
            scalar a, b, c, d;
            scalar Fcent = ((1-a) exp(-T/b) + a exp (-T/c) + exp (-d/T));
        }*/
    }

    //- Move calculated reaction rates into chemData_::reacRates_
    chemData_.update_k(k);
}


AFC::scalar AFC::Chemistry::calcM
(
    const map<word, scalar>& speciesMol,
    const int& r
)
{
    const wordList& species = chemData_.species();

    const wordList& Mcomp = chemData_.Mcomp(r);
    forAll (Mcomp, s)
    {
        Info << r << ": species: " << s << " -> " << Mcomp.at(s) << endl;
    }

/*    forAll(speciesMol, i)
    {
        Info<< "speciesMol:" << species[i] << " -- " <<speciesMol.at(species[i]) << endl;
    }*/
    return 1;
}


AFC::scalar AFC::Chemistry::calcM_Warnatz
(
    const map<word, scalar>& C
)
{
    return
    (
        C.at("H2")
      + 6.5*C.at("H2O")
      + 0.4*C.at("O2")
      + 0.4*C.at("N2")
      + 0.75*C.at("CO")
      + 1.5*C.at("CO2")
      + 3.0*C.at("CH4")
    );
}


AFC::scalar AFC::Chemistry::arrhenius
(
    const scalar& A,
    const scalar& beta,
    const scalar& Ea,
    const scalar& T
)
{
    return (A * pow(T, beta) * exp(Ea/(AFC::Constants::R*T)));
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
