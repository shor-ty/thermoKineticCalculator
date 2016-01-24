/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2015 Holzmann-cfd
  c     c-o-o-o             |
------------------------------------------------------------------------------
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

#include "chemistryCalc.hpp"
#include "constants.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::ChemistryCalc::ChemistryCalc()
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::ChemistryCalc::~ChemistryCalc()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * * //

void AFC::ChemistryCalc::updatekfkb
(
    const scalar& T,
    const Thermo& thermo,
    ChemistryData& chemData
)
{
    //- No. of reactions
    const int& reac = chemData.nReac();

    //- Calculate reaction rates k
    for (int r=0; r<=reac; r++)
    {
        //- Take backward reaction into account?
        bool BW{false};

        //- Backward reaction included
        if (chemData.BR(r))
        {
            BW = true; 
        }

        //- a) Update kf
        updateKf(r, T, chemData);

        //- b) if reversible reaction, backward reaction is needed
        //  therefore we need free gibbs energy
        if (BW)
        {
            //- 1) Update Kc
            updateKc(r, T, thermo, chemData);

            //- 2) Update kb
            //updateKb();
        }
    }
}


void AFC::ChemistryCalc::updateKf
(
    const int& r,
    const scalar& T,
    ChemistryData& chemData
)
{
    //- Arrhenius coeffs
    const scalarField& arrCoeffs = chemData.arrheniusCoeffs(r);

    //- Pre-exponential factor, Unit depend on reaction
    //  + unimolecular reaction [1/s]
    //  + bimolecular reaction [cm^3/mol/s]
    //  + trimolecular reaction [cm^6/mol^2/s]
    const scalar A = arrCoeffs[0];

    //- Expontent
    const scalar beta = arrCoeffs[1];

    //- Activation energy [cal/mol]
    const scalar Ea  = arrCoeffs[2];  

    //- Stardard elementar reaction without any TBR
    //  UNITS:
    //  + uni-moleculare reaction: [1/s]
    //  + bi--moleculare reaction: [cm^3/mol/s]
    //  + tri-moleculare reaction: [cm^6/mol^2/s]
    chemData.updateKf
    (
        r,
        this->arrhenius(A, beta, Ea, T)
    );
}


void AFC::ChemistryCalc::updateKb
(
    const int& r,
    ChemistryData& chemData
)
{
    //- Update kb
    chemData.updateKb
    (
        r,
        chemData.kf(r) / chemData.Kc(r)
    );
}


void AFC::ChemistryCalc::updateKc
(
    const int& r,
    const scalar& T,
    const Thermo& thermo,
    ChemistryData& chemData
)
{
    //- Calculate sum of entropy of species in reaction r
    scalar deltaH{0};
    scalar deltaS{0};

    //- Scalar list of stochiometric factors of reaction r
    const scalarList& nu = chemData.nu(r);

    {
        const wordList& species = chemData.speciesInReaction(r); 

        forAll(species, s)
        {
            deltaH += thermo.H(species[s], T) * nu[s];
            deltaS += thermo.S(species[s], T) * nu[s];
        }
    }

    const scalar deltaG = deltaH - deltaS * T ;

    //- Calculate reaction rate constant Kp
    const scalar Kp = exp(-1 * deltaG / (AFC::Constants::Rcal * T ));

    //- Calculate exponent
    int exponent{0}; 

    forAll(nu, i)
    {
        exponent += nu[i];
    }

    //- Get the pressure of the calculation
    const scalar& p = thermo.p();

    //- Calculate the equilibrium constant Kc
    //  + Pressure in atm (maybe always 1e6?)
    const scalar Kc = Kp * pow(exp(p*10/AFC::Constants::Rcal/T), exponent);

    //- Update Kc
    chemData.updateKc(r, Kc);
}
    

AFC::scalar AFC::ChemistryCalc::calcM_Warnatz
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


AFC::scalar AFC::ChemistryCalc::arrhenius
(
    const scalar& A,
    const scalar& beta,
    const scalar& Ea,
    const scalar& T
) const
{
    return
    (
        A * pow(T, beta) *
        exp( Ea / (AFC::Constants::Rcal * T) ) 
    );
}


void AFC::ChemistryCalc::omega
(
    const map<word, scalar>& speciesMol,
    ChemistryData& chemData
)
{
    //- Species list
    const wordList& species = chemData.species();

    //- Calculate omega for each species
    forAll(species, s)
    {
        //- Reaction no. where the species[s] is invoken
        //const scalarField& reacNo = chemData.reacNoForSpecies(s);
        scalar t =  speciesMol.at("H");

        t = t+1;
    }
}


// ************************************************************************* //
