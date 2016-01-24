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
    const map<word, scalar>& speciesCon,
    const Thermo& thermo,
    ChemistryData& chemData
)
{
    //- No. of reactions
    const int& reac = chemData.nReac();

    //- Calculate reaction rates kf and kb (<- if needed)
    for (int r=0; r<reac; r++)
    {
        //- a) Check if third body reaction and if yes, which one
        //  + enhanced -> false, calculated [M] with all species
        //  + enhanced -> true, calculated [M] with enhanced factors
        bool enhanced = false;

        const bool TBR = thirdBodyReaction(r, enhanced, chemData);

        scalar M = 0;

        //- Calculate [M] if needed
        if (TBR)
        {
            M = calculateM(r, speciesCon, chemData);
        }

        //- Take backward reaction into account?
        bool BW{false};

        //- Backward reaction included
        if (chemData.BR(r))
        {
            BW = true; 
        }

        //- a) Update kf
        updateKf(r, T, M, chemData);

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
    const scalar& M,
    ChemistryData& chemData
)
{
    //- Reaction rate kf
    scalar kf = 0;

    //- Arrhenius coeffs
    const scalarField& arrCoeffs = chemData.arrheniusCoeffs(r);

    //- Pre-exponential factor, Unit depend on reaction
    //  + unimolecular reaction [1/s]
    //  + bimolecular reaction [cm^3/mol/s]
    //  + trimolecular reaction [cm^6/mol^2/s]
    const scalar& A = arrCoeffs[0];

    //- Expontent
    const scalar& beta = arrCoeffs[1];

    //- Activation energy [cal/mol]
    const scalar& Ea  = arrCoeffs[2];  

    //- Lindemann (LOW) | TROE or SRI
    if (chemData.LOW(r))
    {
        //- a) calculate kinf with normal (high) arrhenius coeffs 
        scalar kinf = arrhenius(A, beta, Ea, T);

        //- b) get arrhenius coeffs for low pressure
        const scalarField& arrCoeffsLow = chemData.LOWCoeffs(r);

            //- Pre-exponental factor
            const scalar& ALow = arrCoeffsLow[0];
            
            //- Exponent
            const scalar& betaLow = arrCoeffsLow[1];

            // Activation energy [cal/mol]
            const scalar& EaLow = arrCoeffsLow[2];

        //- c) calculate k0 with low coeffs 
        scalar k0 = arrhenius(ALow, betaLow, EaLow, T);

        //- d) calculate reduced pressure pred
        scalar pred = k0 * M / kinf; 

        //- e) Calculate F
        //  + Complex functions if TROE or SRI
        //  + For Lindemann = 1
        scalar F = 1;

        //- TROE formula for F
        if (chemData.TROE(r))
        {
            //- i) Get TROE coeffs
            const scalarField& troeCoeffs = chemData.TROECoeffs(r);

            //- Alpha
            const scalar& alpha = troeCoeffs[0];

            //- T***
            const scalar& Tsss = troeCoeffs[1];

            //- T**
            const scalar& Tss = troeCoeffs[2];

            //- T*
            const scalar& Ts = troeCoeffs[3];
            
            //- ii) calculate Fcent
            const scalar Fcent = (1 - alpha) * exp(-1 * T / Tsss) +
                alpha * exp(-1 * T / Ts) + exp(-1 * Tss / T);
            
            //- iii) calculate constants
            const scalar c = -0.4 - 0.67 * log(Fcent);
            const scalar n = 0.75 - 1.27 * log(Fcent);
            const scalar d = 0.14;

            //- iv) calculate logF
            const scalar logF = 
                pow((1 + pow((log(pred) + c)/(n - d*(log(pred)+c)),2)), -1) * log(Fcent);

            //- v) get F
            F = exp(logF);
        }
        //- SRI formula for F
        if (chemData.SRI(r))
        {
            //- i) Get SRI coeffs
            const scalarField& sriCoeffs = chemData.SRICoeffs(r);

            //- a
            const scalar a = sriCoeffs[0];

            //- b
            const scalar b = sriCoeffs[1];

            //- c
            const scalar c = sriCoeffs[2];

            //- d
            const scalar d = sriCoeffs[3];

            //- e
            const scalar e = sriCoeffs[4];

            //- ii) calculate X
            const scalar X = 1 / (1 + pow(log(pred), 2));

            //- iii) calculate F
            F = d * pow(a * exp(-1 * b / T) + exp (-1 * T / c), X) * pow(T, e);
        } 

        //- f) calculate kf
        kf = kinf * (pred / (1 + pred)) * F;
    }
    //- Normal reaction
    else
    {
        kf = arrhenius(A, beta, Ea, T);
    }

    //- Update
    //  UNITS:
    //  + uni-moleculare reaction: [1/s]
    //  + bi--moleculare reaction: [cm^3/mol/s]
    //  + tri-moleculare reaction: [cm^6/mol^2/s]
    chemData.updateKf
    (
        r,
        kf
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


AFC::scalar AFC::ChemistryCalc::calculateM
(
    const int& r,
    const map<word, scalar>& speciesCon,
    const ChemistryData& data 
)
{
    //- Tmp M
    scalar M{0};

    //- ENHANCED species
    wordList enhancedSpecies;

    const bool& EH = data.ENHANCED(r);


    //- Use enhance values
    if (EH)
    {
        //- Get information
        enhancedSpecies = data.enhancedSpecies(r);
    }

    //- Each species is used as third body partner
    {
        const wordList& species = data.species();

        //- [M] = sum of all concentrations
        forAll(species, s)
        {
            if (EH)
            {
                //- Search if Species available
                forAll(enhancedSpecies, e)
                {
                    //- If found, then use modified value
                    if (enhancedSpecies[e] == species[s])
                    {
                        M += speciesCon.at(species[s]) * data.enhancedFactors(r, species[s]);
                    }
                    //- otherwise use value 1
                    else
                    {
                        M += speciesCon.at(species[s]);
                    }
                }
            }
            else
            {
                M += speciesCon.at(species[s]);
            }
        }
    }

    //- Return M [mol/m^3]
    return M;
}


bool AFC::ChemistryCalc::thirdBodyReaction
(
    const int& r,    
    bool& enhanced,
    const ChemistryData& data
)
{
    //- Get third body information
    const bool& TBR = data.TBR(r);

    //- Get information about enhanced factors
    enhanced = data.ENHANCED(r);

    return TBR;
}


// ************************************************************************* //
