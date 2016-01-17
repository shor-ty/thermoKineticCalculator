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

void AFC::ChemistryCalc::k
(
    const scalar& T,
    const map<word,scalar>& speciesMol,
    const Thermo& thermo,
    ChemistryData& chemData
)
{
    //- No. of reactions
    const int& reac = chemData.nReac();

    //- Store forward reaction rates kf (tmp)
    scalarField kf;

    //- Store forward reaction rates kb (tmp)
    scalarField kb;

    //- REMOVE TODO
    scalar t =  speciesMol.at("H");
    t = t+1;

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

        //- a) Calculate kf
        calculateKf(kf, r, T, speciesMol, chemData);

        //- b) if reversible reaction, backward reaction is needed
        //  therefore we need free gibbs energy
        if (BW)
        {
            calculateKb(kb, kf, r, T, speciesMol, thermo, chemData);
        }
    }

    //- Move calculated forward reaction rates into chemData_::kf_
    chemData.updateKf(kf);
}


void AFC::ChemistryCalc::calculateKf
(
    scalarField& kf,
    const int& r,
    const scalar& T,
    const map<word, scalar>& speciesMol,
    ChemistryData& chemData
)
{

    //- Arrhenius coeffs
    const scalarField& arrCoeffs = chemData.arrheniusCoeffs(r);

    //- Pre-exponential factor [cm^3/mol/s]
    //  combustion.berkeley.edu
    const scalar A = arrCoeffs[0];

    //- Expontent
    const scalar beta = arrCoeffs[1];

    //- Activation energy [cal/mol]
    const scalar Ea  = arrCoeffs[2];  

    //- Stardard elementar reaction without any TBR
    kf.push_back(this->arrhenius(A, beta, Ea, T));
}


void AFC::ChemistryCalc::calculateKb
(
    scalarField& kb,
    const scalarField& kf,
    const int& r,
    const scalar& T,
    const map<word, scalar>& speciesMol,
    const Thermo& thermo,
    ChemistryData& chemData
)
{
    //- Calculate free Gibbs energy for reaction r
//    const scalar G = thermo.G(r, T);
    scalar F = 1;

    kb.push_back
    (
        kf[r] / exp(F / (AFC::Constants::R * T ))
    );
}

AFC::scalar AFC::ChemistryCalc::calcM
(
    const map<word, scalar>& speciesMol,
    const int& r
)
{
//    const wordList& species = this->species();

//    const wordList& Mcomp = chemData_.Mcomp(r);

/*    forAll(speciesMol, i)
    {
        Info<< "speciesMol:" << species[i] << " -- " <<speciesMol.at(species[i]) << endl;
    }*/
    // TODO
    return speciesMol.at("H2") * r;
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
)
{
    return (A * pow(T, beta) * exp(Ea/(AFC::Constants::R*T)));
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
