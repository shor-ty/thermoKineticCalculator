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
    if (debug)
    {
        Info<< "Constructor ChemistryCalc\n" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::ChemistryCalc::~ChemistryCalc()
{
    if (debug)
    {
        Info<< "Destructor ChemistryCalc\n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * * //

AFC::scalar AFC::ChemistryCalc::kf
(
    const int& r,
    const scalar& T,
    const ChemistryData& chemData,
    const bool LOW
) const
{
    //- Pre-exponential factor, Unit depend on reaction
    //  + unimolecular reaction [1/s]
    //  + bimolecular reaction [cm^3/mol/s]
    //  + trimolecular reaction [cm^6/mol^2/s]

    //- TODO make faster with pointer or reference
    scalarField arrheniusCoeffs;
   
    if (LOW)
    {
        arrheniusCoeffs = chemData.LOWCoeffs(r);
    } 
    else
    {
        arrheniusCoeffs = chemData.arrheniusCoeffs(r);
    }

    return arrhenius
        (
            arrheniusCoeffs[0],
            arrheniusCoeffs[1],
            arrheniusCoeffs[2],
            T
        );
}


AFC::scalar AFC::ChemistryCalc::kb
(
    const int& r,
    const scalar& T,
    const ChemistryData& chemData,
    const Thermo& thermo,
    const bool LOW
) const
{
    //- Calculate backward reaction rate kb
    return (kf(r, T, chemData, LOW) / keq(r, T, chemData, thermo));
}


AFC::scalar AFC::ChemistryCalc::keq
(
    const int& r,
    const scalar& T,
    const ChemistryData& chemData,
    const Thermo& thermo
) const
{
    //- Calculate sum of free GIBBS energy of reaction r
    const scalar deltaG = dG(r, T, chemData, thermo);

    //- Calculate reaction rate constant kp
    const scalar Kp = exp(-1 * deltaG / (AFC::Constants::R * T ));

    //- Global reaction order == exponent
    scalar exponent = chemData.globalReactionOrder(r);

    //- Save some calculation because if exponent == 0 hence keq == kp
    if (exponent == 0)
    {
        return Kp;
    }
    else
    {
        //- Get the pressure of the calculation
        const scalar& p = thermo.p();


        //- Calculate the equilibrium constant Keq and return it
        //  In the formulation p is normally in [dyn/cm^2] and R in
        //  [ergs/K/mol]. We use [Pa] and [J/K/mol] therefore the
        //  have to multiply R by 1000 | p/R for 1 bar ~ 12.1802
        //  TODO change R to Rcal
        return (Kp / pow(AFC::Constants::Rcal * 1e3 / p * T, exponent));
    }
}
    

AFC::scalar AFC::ChemistryCalc::arrhenius
(
    const scalar& A,
    const scalar& beta,
    const scalar& Ea,
    const scalar& T
) const
{
    //-TODO change R to RCalc again R -> must be [cal]
    return
    (
        A * pow(T, beta) * exp(-1*Ea / (AFC::Constants::R * T))
    );
}


AFC::scalar AFC::ChemistryCalc::Fcent
(
    const int& r,
    const scalar& T,
    const ChemistryData& chemData
) const
{
    //- TODO speedup
    //- TROE coefficients
    const scalarField& troeCoeffs = chemData.TROECoeffs(r);

    //- Alpha
    const scalar& alpha = troeCoeffs[0];

    //- T***
    const scalar& Tsss = troeCoeffs[1];

    //- T*
    const scalar& Ts = troeCoeffs[2];

    //- T**
    const scalar& Tss = troeCoeffs[3];
    
    //- Return Fcent 
    return ((1 - alpha) * exp(-1 * T / Tsss) +
        alpha * exp(-1 * T / Ts) + exp(-1 * Tss / T));
}


AFC::scalar AFC::ChemistryCalc::Flog
(
    const int& r,
    const scalar& T,
    const scalar& M,
    const ChemistryData& chemData
) const
{
    //- a) Calculate Fcent
    const scalar& F_cent = Fcent(r, T, chemData);

    //- b) Calculate constants
    const scalar N = 0.75 - 1.27 * log(F_cent);
    const scalar C = -0.4 - 0.67 * log(F_cent);

    //- c) Calculate reduced pressure Pr
    //  + Here we need k for LOW and normal pressure
    const scalar& kinf = kf(r, T, chemData);
    const scalar& klow = kf(r, T, chemData, true);

    const scalar Pr = klow * M / kinf;

    //- d) Calculate enumerator
    const scalar& enumerator = log10(F_cent);

    //- e) Calculate denominator
    const scalar& denominator = 1 +
        pow((log10(Pr) + C)/(N - 0.14 * (log10(Pr) + C)),2);

    //- f) Return Flog
    return (enumerator/denominator);
}


AFC::scalar AFC::ChemistryCalc::dH
(
    const int& r,
    const scalar& T,
    const ChemistryData& data,
    const Thermo& thermo
) const
{
    //- Stochiometric factors of reaction r
    const map<word, scalar>& educts = data.nuEducts(r);
    const map<word, scalar>& products = data.nuProducts(r);

    scalar dH{0};

    forAll(educts, e)
    {
        dH += thermo.H(e.first, T) * e.second;
    }

    forAll(products, p)
    {
        dH += thermo.H(p.first, T) * p.second;
    }

    return dH;
}


AFC::scalar AFC::ChemistryCalc::dG
(
    const int& r,
    const scalar& T,
    const ChemistryData& data,
    const Thermo& thermo
) const
{
    //- Stochiometric factors of reaction r
    const map<word, scalar>& educts = data.nuEducts(r);
    const map<word, scalar>& products = data.nuProducts(r);

    scalar dG{0};

    forAll(educts, e)
    {
        dG += thermo.G(e.first, T) * e.second;
    }

    forAll(products, p)
    {
        dG += thermo.G(p.first, T) * p.second;
    }

    return dG;
}


AFC::scalar AFC::ChemistryCalc::dS
(
    const int& r,
    const scalar& T,
    const ChemistryData& data,
    const Thermo& thermo
) const
{
    //- Stochiometric factors of reaction r
    const map<word, scalar>& educts = data.nuEducts(r);
    const map<word, scalar>& products = data.nuProducts(r);

    scalar dS{0};

    forAll(educts, e)
    {
        dS += thermo.S(e.first, T) * e.second;
    }

    forAll(products, p)
    {
        dS += thermo.S(p.first, T) * p.second;
    }

    return dS;
}


// ************************************************************************* //

    /*if (chemData.LOW(r))
    {
        //- a) calculate kinf with normal (high) arrhenius coeffs 
        const scalar& kinf = arrhenius(A, beta, Ea, T);

        //- b) get arrhenius coeffs for low pressure
        const scalarField& arrCoeffsLow = chemData.LOWCoeffs(r);

        //- c) calculate k0 with low coeffs 
        const scalar k0 =
            arrhenius
            (
                arrCoeffsLow[0],
                arrCoeffsLow[1],
                arrCoeffsLow[2],
                T
            );

        //- d) calculate reduced pressure pReduced
        const scalar pReduced = k0 * chemData.M() / kinf; 

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
                pow((1 + pow((log(pred) + c)/(n - 0.14*(log10(pred)+c)),2)), -1) * log10(Fcent);

            //- v) get F
            F = pow(10, logF);

            if (r == 12)
            {
                Info<< "Fcent = (1 - " << alpha << ") * e^(-" << T << "/"
                 << Tsss << ") + " << alpha << " * e^(-" << T << "/" << Ts
              << ") + e^(-" << Tss << "/" << T << ")\n";   
                Info<< "Fcent: " << Fcent << endl;
                Info<< "N: " << n << "\n";
                Info<< "C: " << c << "\n";
                Info<< "k0: " << k0 << "\n";
                Info<< "kinf: " << kinf << "\n";
                Info<< "[M]: "<< chemData.M() << "\n";
                Info<< "Pred: " << pred << "\n";
                Info<< "logF: " << logF << "\n";
                Info<< "F: " << F << "\n";
            }
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
        chemData.updateKf(kinf * (pred / (1 + pred)) * F);
    }
    
    //- Normal reaction
    //else
    {
        return (arrhenius(A, beta, Ea, T));
    }

    //- Return
    //  UNITS:
    //  + uni-moleculare reaction: [1/s]
    //  + bi--moleculare reaction: [cm^3/mol/s]
    //  + tri-moleculare reaction: [cm^6/mol^2/s]

}*/
/*void AFC::ChemistryCalc::calculatekfkb
(
    const int& r,
    const scalar& T,
    const map<word, scalar>& speciesCon,
    const Thermo& thermo,
    ChemistryData& chemData
)
{
    const bool TBR = thirdBodyReaction(r, chemData);

    //- Calculate [M] if needed and store
    if (TBR)
    {
        calculateM(r, speciesCon, chemData);
    }

    //- Set third body reaction concentration [M] [g/cm^3]
    //chemData.setM(M);

    //- Take backward reaction into account?
    //const bool& BW = chemData.BR(r);

    //- a) calculate kf
    calculateKf(r, T, chemData);

    //- b) if reversible reaction, backward reaction is needed
    //  therefore we need free gibbs energy
    
    //if (BW)
    {
        //- 1) calculate Kc
        const scalar Kc = calculateKc(r, T, thermo, chemData);

        //- 2) calculate kb
        kb = calculateKb(kf, Kc);
    }
}


AFC::scalar AFC::ChemistryCalc::kf
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
        const scalar pred = k0 * chemData.M() / kinf; 

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
                pow((1 + pow((log(pred) + c)/(n - 0.14*(log10(pred)+c)),2)), -1) * log10(Fcent);

            //- v) get F
            F = pow(10, logF);

            if (r == 12)
            {
                Info<< "Fcent = (1 - " << alpha << ") * e^(-" << T << "/"
                 << Tsss << ") + " << alpha << " * e^(-" << T << "/" << Ts
              << ") + e^(-" << Tss << "/" << T << ")\n";   
                Info<< "Fcent: " << Fcent << endl;
                Info<< "N: " << n << "\n";
                Info<< "C: " << c << "\n";
                Info<< "k0: " << k0 << "\n";
                Info<< "kinf: " << kinf << "\n";
                Info<< "[M]: "<< chemData.M() << "\n";
                Info<< "Pred: " << pred << "\n";
                Info<< "logF: " << logF << "\n";
                Info<< "F: " << F << "\n";
            }
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
        chemData.updateKf(kinf * (pred / (1 + pred)) * F);
    }
    //- Normal reaction
    else
    {
        chemData.updateKf(arrhenius(A, beta, Ea, T));
    }

    //- Return
    //  UNITS:
    //  + uni-moleculare reaction: [1/s]
    //  + bi--moleculare reaction: [cm^3/mol/s]
    //  + tri-moleculare reaction: [cm^6/mol^2/s]
//}
*/

/*
AFC::scalar AFC::ChemistryCalc::calculateOmega
(
    const word& species,
    const scalar& T,
    map<word, scalar>& con,
    const Thermo& thermo,
    ChemistryData& chemData
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

        //- Elementarreaction no.
        const unsigned int& r = inReaction[i];

        //- Temporar fields
        scalar kf{0};
        scalar kb{0};

        //- Calculate reaction rate kf and kb for reaction r
        calculatekfkb(r, T, con, thermo, chemData);

        //- Get all species that are used in reaction r
        const wordList& speciesInReaction = chemData.speciesInReaction(r);


        //- Get all species that are acting as products in reaction r
        const wordList& prodSpecies = chemData.prodSpecies(r);

        //- Get all species that act as educts in reaction r
        const wordList& educSpecies = chemData.educSpecies(r);

        //- Scalar list of stochiometric factors of reaction r
        const scalarList& nu = chemData.nu(r);

        //- Stochiometric factors of products
        map<word, scalar> nuProd = chemData.nuProd(r);

        //- Stochiometric factors of educts
        map<word, scalar> nuEduc = chemData.nuEduc(r);

        //- TMP fields
        scalar prod{1};
        scalar educ{1};


        //- Product side, con is in [mol/cm^3]
        forAll(prodSpecies, s)
        {
            prod *= pow(con.at(prodSpecies[s]), nuProd.at(prodSpecies[s]));
        }

        //- Educt side, con is in [mol/cm^3]
        forAll(educSpecies, s)
        {
            educ *= pow(con.at(educSpecies[s]), nuEduc.at(educSpecies[s]));
        }

        const scalar nuSpecies = nuProd[species] - nuEduc[species];

        //- kf is in cm^3, prod and educ in mol/cm^3
        const scalar& kf = chemData.kf();
        const scalar& kb = chemData.kb();

        omega += chemData.M() * nuSpecies * (kf * educ + kb * prod);
    } 

    return omega;
}


void AFC::ChemistryCalc::calculateM
(
    const int& r,
    const map<word, scalar>& speciesCon,
    ChemistryData& data 
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

    //- Store M [mol/cm^3]
    chemData.updateM(M);
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
*/

