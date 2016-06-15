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

void AFC::ChemistryCalc::calculatekfkb
(
    const int& r,
    scalar& kf,
    scalar& kb,
    const scalar& T,
    const map<word, scalar>& speciesCon,
    const Thermo& thermo,
    const ChemistryData& chemData
)
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
    const bool& BW = chemData.BR(r);

    //- a) calculate kf
    kf = calculateKf(r, T, M, chemData);

    //- b) if reversible reaction, backward reaction is needed
    //  therefore we need free gibbs energy
    
    if (BW)
    {
        //- 1) calculate Kc
        const scalar Kc = calculateKc(r, T, thermo, chemData);

        //- 2) calculate kb
        kb = calculateKb(kf, Kc);
    }
}


AFC::scalar AFC::ChemistryCalc::calculateKf
(
    const int& r,
    const scalar& T,
    const scalar& M,
    const ChemistryData& chemData
)
{
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
        const scalar k0 = arrhenius(ALow, betaLow, EaLow, T);

        //- d) calculate reduced pressure pred
        const scalar pred = k0 * M / kinf; 

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

            //- T*
            const scalar& Ts = troeCoeffs[2];

            //- T**
            const scalar& Tss = troeCoeffs[3];
            
            //- ii) calculate Fcent
            const scalar Fcent = (1 - alpha) * exp(-1 * T / Tsss) +
                alpha * exp(-1 * T / Ts) + exp(-1 * Tss / T);
            
            //- iii) calculate constants
            const scalar c = -0.4 - 0.67 * log(Fcent);
            const scalar n = 0.75 - 1.27 * log(Fcent);
            const scalar d = 0.14;

            //- iv) calculate logF
            const scalar logF = 
                pow((1 + pow((log(pred) + c)/(n - d*(log10(pred)+c)),2)), -1) * log10(Fcent);

            //- v) get F
            F = pow(10, logF);
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
        return kinf * (pred / (1 + pred)) * F;
    }
    //- Normal reaction
    else
    {
        return arrhenius(A, beta, Ea, T);
    }

    //- Return
    //  UNITS:
    //  + uni-moleculare reaction: [1/s]
    //  + bi--moleculare reaction: [cm^3/mol/s]
    //  + tri-moleculare reaction: [cm^6/mol^2/s]
}


AFC::scalar AFC::ChemistryCalc::calculateKb
(
    const scalar& kf,
    const scalar& Kc
)
{
    //- Update kb
    return (kf / Kc);
}


AFC::scalar AFC::ChemistryCalc::calculateKc
(
    const int& r,
    const scalar& T,
    const Thermo& thermo,
    const ChemistryData& chemData
)
{
    //- Calculate sum of entropy of species in reaction r
    scalar deltaH{0};
    scalar deltaS{0};

    //- Scalar list of stochiometric factors of reaction r
    const map<word, scalar>& nu = chemData.nu(r);

    {
        //const wordList& species = chemData.speciesInReaction(r); 

        /*forAll(nu, s)
        {
            deltaH += thermo.H(s, T) * nu.at(s)
            deltaS += thermo.S(s, T) * nu.at(s);
        }*/
    }

    const scalar deltaG = deltaH - deltaS * T ;

    //- Calculate reaction rate constant Kp
    const scalar Kp = exp(-1 * deltaG / (AFC::Constants::R * T ));

    //- Calculate exponent
    int exponent{0}; 

    forAll(nu, s)
    {
//        exponent += nu.at(i);
    }

    //- Get the pressure of the calculation
    const scalar& p = thermo.p();

    //- Calculate the equilibrium constant Kc and return it
    //  In the formulation p is normally in [dyn/cm^2] and R in [ergs/K/mol]
    //  We use [Pa] and [J/K/mol] therefore the have to multiply p * 10 and
    //  R by 10e7
    scalar Kc = (Kp * pow((p*10/AFC::Constants::R/1e7/T), exponent));

    return Kc;
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


void AFC::ChemistryCalc::calculateOmega
(
    const scalar& T,
    map<word, scalar>& con,
    const Thermo& thermo,
    const ChemistryData& chemData
)
{
    //- Species list
    const wordList& species = chemData.species();

    //- Calculate Omega using GAUSS-SEIDEL
    //  Each species has to be calculated
    forAll(species, s)
    {
        //- Get reaction no. where species is included
        //const List<int>& inReaction = chemData.reacNumbers(species[s]); 

        //forAll(inReaction, r)
        {
            //- Temporar fields
            scalar kf{0};
            scalar kb{0};

            //- Calculate reaction rate kf and kb for reaction i
        //    calculatekfkb(r, kf, kb, T, con, thermo, chemData);
        }
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
                //- Species in enhanced found?
                bool found{false};

                //- Search if Species available
                forAll(enhancedSpecies, e)
                {
                    //- If found, then use modified value
                  //  if (enhancedSpecies[e] == species[s])
                    {
                        //- Species found, use modified value and skip using value 1
                        found = true;
                        //M += speciesCon.at(species[s]) * data.enhancedFactors(r, species[s]);
                    }
                }

                if (!found)
                {
                    //M += speciesCon.at(species[s]);
                }
            }
            else
            {
                //M += speciesCon.at(species[s]);
            }
        }
    }

    //- Return M [mol/cm^3]
    //  M is in mol/m^3 therefore divide by 1000000 cm^3/m^3
    return (M / 1000000);
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
