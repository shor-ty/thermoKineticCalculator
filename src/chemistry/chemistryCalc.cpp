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
    const ChemistryData& chemData
) const
{
    int r = 1;
    
    //- Arrhenius coeffs
    const scalarField& arrCoeffs = chemData.arrheniusCoeffs(r);

    //- Pre-exponential factor [cm^3/mol/s]
    //  combustion.berkeley.edu
    const scalar A = arrCoeffs[0];

    //- Expontent
    const scalar beta = arrCoeffs[1];

    //- Activation energy [cal/mol]
    const scalar Ea  = arrCoeffs[2];  

    const scalar T{300};

    scalar kf = this->arrhenius(A, beta, Ea, T);

    Info<< "kf: " << kf << endl;

    scalar H = -104200; //- cal

    scalar S = -23.6; //- cal / mol K

    scalar G = H - S * T;

    scalar p = 1e6;

    scalar Kp = exp(-1 * G / (AFC::Constants::Rcal * T));

    scalar Kc = Kp * pow((p/(8.3144598e7 * T)), -1);

    Info<< "Kp = " << Kp << " --> K: " << (kf / Kp) << endl;
    Info<< "Kc = " << Kc << " --> K: " << (kf / Kc) << endl;

}

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
    scalar G;

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
    //  For SI units cal -> joule
    chemData.updateKf
    (
        r,
        this->arrhenius(A, beta, Ea, T) * AFC::Constants::calToJoule
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
    //- Get free Gibbs energy G
    //const scalar& G = thermo.G();

    scalar G{0};

    //- Calculate reaction rate constant Kp
    const scalar Kp = exp(G / (AFC::Constants::R * T ));
        
    //- Get stochiometric factors
    const scalarList& nu = chemData.nu(r);

    //- Calculate exponent
    int exponent{0}; 

    forAll(nu, i)
    {
        exponent += nu[i];
    }

    scalar p{0};

    const scalar Kc = Kp * pow(exp(p/AFC::Constants::R/T), exponent);

    //- Update Kc
    chemData.updateKc(r, Kc);
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
) const
{
    //- Unit depend on reaction
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
