/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     T hermo       | Open Source Thermo-Kinetic Library
  c-o-o-c     K iknetic     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
  c     c-o-o-o             |
------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Constructor.

    TKC is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    TKC is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TKC; if not, see <http://www.gnu.org/licenses/>

\*---------------------------------------------------------------------------*/

#include "chemistryCalc.hpp"
#include "constants.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::ChemistryCalc::ChemistryCalc(const string fileName, const Thermo& thermo)
:
    ChemistryData(fileName),
    thermo_(thermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::ChemistryCalc::~ChemistryCalc()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const TKC::Thermo TKC::ChemistryCalc::thermo() const
{
    return thermo_;
}


// * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * * //

TKC::scalar TKC::ChemistryCalc::kf
(
    const int r,
    const scalar T,
    const map<word, scalar>& c,
    const bool kb
) const
{
    //- Forward reaction possible?
    //  If we come from kb, also calculate kf
    if (forwardReaction(r) || kb)
    {
        //- Pre-exponential factor, Unit depend on reaction
        //  + unimolecular reaction [1/s]
        //  + bimolecular reaction [cm^3/mol/s]
        //  + trimolecular reaction [cm^6/mol^2/s]
        const scalarField& arrCoeffs = arrheniusCoeffs(r);

        //- Is a third body collision partner included?
        if (TBR(r))
        {
            //- If no Lindemann formula is used, standard arrhenius
            if (!LOW(r) && !TROE(r) && !SRI(r))
            {
                //- Calculate standard arrhenius
                return arrhenius(arrCoeffs[0], arrCoeffs[1], arrCoeffs[3], T);
            }
            //- Lindemann formulation
            else if (LOW(r) && !TROE(r))
            {
                return Lindemann(r, T, c);
            }
            //- TROE formulation
            else if (LOW(r) && TROE(r) && !SRI(r))
            {

            }
            //- SRI formulation
            else if (LOW(r) && TROE(r) && SRI(r))
            {

            }
            else
            {
                Warning
                (
                    "Reaction " + elementarReaction(r) + " is not\n    "
                    "defined regarding the fall off regions. No option fits\n",
                    __FILE__,
                    __LINE__
                );
            }
        }
        //- Standard reaction
        else
        {
            return arrhenius(arrCoeffs[0], arrCoeffs[1], arrCoeffs[2], T);
        }
    }
    else
    {
        return 0;
    }

    //- Compress compiler warning
    return 0;
}


TKC::scalar TKC::ChemistryCalc::kb
(
    const int r,
    const scalar T,
    const map<word, scalar>& c
) const
{
    //- Backward reaction possible?
    if (backwardReaction(r))
    {
        //- Calculate backward reaction rate kb
        return (kf(r, T, c, true) / keq(r, T));
    }
    else
    {
        return 0;
    }
}


TKC::scalar TKC::ChemistryCalc::keq
(
    const int r,
    const scalar T
) const
{
    //- Calculate sum of free GIBBS energy of reaction r
    const scalar deltaG = dg(r, T);

    //- Calculate reaction rate constant kp
    const scalar Kp = exp(-1 * deltaG / (TKC::Constants::R * T ));

    //- Global reaction order == exponent
    const scalar exponent = globalReactionOrder(r);

    //- Save some calculation because if exponent == 0 hence keq == kp
    if (exponent == 0)
    {
        return Kp;
    }
    else
    {
        //- Get the pressure of the calculation
        const scalar p = thermo_.p();

        //- Calculate the equilibrium constant Keq and return it
        //  In the formulation p is normally in [dyn/cm^2] and R in
        //  [ergs/K/mol]. We use [Pa] and [J/K/mol] therefore the
        //  have to multiply R by 1000 | p/R for 1 bar ~ 12.1802
        //  TODO change R to Rcal
        return (Kp / pow(TKC::Constants::Rcal * 1e3 / p * T, exponent));
    }
}


TKC::scalar TKC::ChemistryCalc::arrhenius
(
    const scalar A,
    const scalar beta,
    const scalar Ea,
    const scalar T
) const
{
    //-TODO change R to RCalc again R -> must be [cal]
    return
    (
        A * pow(T, beta) * exp(-1*Ea / (TKC::Constants::Rcal * T))
    //    A*pow(T, beta) * exp(-Ea / T)
    );
}


TKC::scalar TKC::ChemistryCalc::Lindemann
(
    const int r,
    const scalar T,
    const map<word, scalar>& c
) const
{
    const scalarField& arrCoeffs = arrheniusCoeffs(r);

    //- Calculate standard arrhenius
    const scalar arrheniusHigh =
        arrhenius(arrCoeffs[0], arrCoeffs[1], arrCoeffs[3], T);

    //- Coefficients for low pressure area
    const scalarField& arrCoeffsLow = LOWCoeffs(r);

    //- Calculate low pressure arrhenius
    const scalar arrheniusLow =
        arrhenius
        (
            arrCoeffsLow[0],
            arrCoeffsLow[1],
            arrCoeffsLow[2],
            T
        );

    //- Calculate the concentration of the third body partner
    const scalar conM = M(r, c);

    //- TODO finish function
    return 0;
}


TKC::scalar TKC::ChemistryCalc::Fcent
(
    const int r,
    const scalar T
) const
{
    //- TROE coefficients
    const scalarField& troeCoeffs = TROECoeffs(r);

    //- Alpha
    const scalar alpha = troeCoeffs[0];

    //- T***
    const scalar Tsss = troeCoeffs[1];

    //- T*
    const scalar Ts = troeCoeffs[2];

    //- T**
    const scalar Tss = troeCoeffs[3];

    //- Return Fcent
    return
    (
        (1 - alpha) * exp(-1 * T / Tsss)
      + alpha * exp(-1 * T / Ts)
      + exp(-1 * Tss / T)
    );
}


TKC::scalar TKC::ChemistryCalc::Flog
(
    const int r,
    const scalar T,
    const scalar M
) const
{
    //- a) Calculate Fcent
    const scalar F_cent = Fcent(r, T);

    //- b) Calculate constants
    const scalar N = 0.75 - 1.27 * log(F_cent);
    const scalar C = -0.4 - 0.67 * log(F_cent);

    //- c) Calculate reduced pressure Pr
    //  + Here we need k for LOW and normal pressure
    const scalar kinf = kf(r, T);
    const scalar klow = kf(r, T);

    const scalar Pr = klow * M / kinf;

    //- d) Calculate enumerator
    const scalar enumerator = log10(F_cent);

    //- e) Calculate denominator
    const scalar denominator = 1 +
        pow((log10(Pr) + C)/(N - 0.14 * (log10(Pr) + C)),2);

    //- f) Return Flog
    return (enumerator/denominator);
}


TKC::scalar TKC::ChemistryCalc::omega
(
    const word species,
    const scalar T,
    const map<word, scalar>& con
) const
{

    //- Calculate Omega for species
    scalar omega{0};

    //- Get reaction no. where species is included
    const List<int>& inReaction = reacNumbers(species);

    //Info<< "Species calculation " << species << "\n"
    //    << "--------------------------------------\n";

    for(size_t i = 0; i<inReaction.size(); ++i)
    {
        //- Elementar Reaction number
        const size_t r = inReaction[i];

        //- Calculate the forward and backward reaction rates
        const scalar kf_ = kf(r, T);
        const scalar kb_ = kb(r, T);

        //- Get all species that are acting as products in reaction r
        const List<word>& prodspecies = products(r);

        //- Get all species that act as educts in reaction r
        const List<word>& educSpecies = educts(r);

        //- Stochiometric factors of products
        map<word, int> nuProd = nuProducts(r);

        //- Stochiometric factors of educts
        map<word, int> nuEduc = nuEducts(r);

        //- TMP fields
        scalar prod{1};
        scalar educ{1};

        //- Product side, con is in [mol/cm^3]
        forAll(prodspecies, s)
        {
            prod *= pow(con.at(s), nuProd.at(s));
        }

        //- Educt side, con is in [mol/cm^3]
        //  Note, abs needed
        forAll(educSpecies, s)
        {
            educ *= pow(con.at(s), abs(nuEduc.at(s)));
        }

        //- Get pre-factor nu'' - nu' :: based on the fact that we already
        //  know the right value, we just have to check if the species
        //  is within the product or educt side
        scalar nuSpecies{0};

        if (nuEduc.count(species) && !nuProd.count(species))
        {
            nuSpecies = nuEduc.at(species);
        }
        else if (!nuEduc.count(species) && nuProd.count(species))
        {
            nuSpecies = nuProd.at(species);
        }
        else if (nuEduc.count(species) && nuProd.count(species))
        {
            nuSpecies = nuProd.at(species) + nuEduc.at(species);
        }
        else
        {
            NotImplemented
            (
                __FILE__,
                __LINE__
            );
        }

        //Info<<std::setw(30) << chemData.elementarReaction(r)
        //    <<std::setw(20) << nuSpecies * (kf_ * educ - kb_ * prod)
        //    << "   "
        //    << nuSpecies << " *( "
        //    << kf_ << " * "
        //    << educ << " - "
        //    << kb_ << " * "
        //    << prod << ")\n" ;
        omega += nuSpecies * (kf_ * educ - kb_ * prod);
    }

    return omega;
}


TKC::scalar TKC::ChemistryCalc::M
(
    const int r,
    const map<word, scalar>& c
) const
{
    //- Tmp M
    scalar M{0};

    const map<word, scalar>& enhancedSpecies = ENHANCEDCoeffs(r);

    //- Each species is used as third body partner
    {
        //- [M] = sum of all concentrations
        loopMap(species, factor, enhancedSpecies)
        {
            M += factor * c.at(species);
        }
    }

    return M;
}

TKC::scalar TKC::ChemistryCalc::dh
(
    const int r,
    const scalar T
) const
{
    //- Stochiometric factors of reaction r
    const map<word, int>& educts = nuEducts(r);
    const map<word, int>& products = nuProducts(r);

    scalar dh{0};

    forAll(educts, e)
    {
        dh += thermo_.h(e.first, T) * e.second;
    }

    forAll(products, p)
    {
        dh += thermo_.h(p.first, T) * p.second;
    }

    return dh;
}


TKC::scalar TKC::ChemistryCalc::dg
(
    const int r,
    const scalar T
) const
{
    //- Stochiometric factors of reaction r
    const map<word, int>& educts = nuEducts(r);
    const map<word, int>& products = nuProducts(r);

    scalar dg{0};

    forAll(educts, e)
    {
        dg += thermo_.g(e.first, T) * e.second;
    }

    forAll(products, p)
    {
        dg += thermo_.g(p.first, T) * p.second;
    }

    return dg;
}


TKC::scalar TKC::ChemistryCalc::ds
(
    const int r,
    const scalar T
) const
{
    //- Stochiometric factors of reaction r
    const map<word, int>& educts = nuEducts(r);
    const map<word, int>& products = nuProducts(r);

    scalar ds{0};

    forAll(educts, e)
    {
        ds += thermo_.s(e.first, T) * e.second;
    }

    forAll(products, p)
    {
        ds += thermo_.s(p.first, T) * p.second;
    }

    return ds;
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
        if (chemData.sRI(r))
        {
            //- i) Get SRI coeffs
            const scalarField& sriCoeffs = chemData.sRICoeffs(r);

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
/*void TKC::ChemistryCalc::calculatekfkb
(
    const int r,
    const scalar& T,
    const map<word, scalar>& speciesCon,
    const Thermo& thermo_,
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
        const scalar Kc = calculateKc(r, T, thermo_, chemData);

        //- 2) calculate kb
        kb = calculateKb(kf, Kc);
    }
}


TKC::scalar TKC::ChemistryCalc::kf
(
    const int r,
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
        if (chemData.sRI(r))
        {
            //- i) Get SRI coeffs
            const scalarField& sriCoeffs = chemData.sRICoeffs(r);

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
void TKC::ChemistryCalc::calculateM
(
    const int r,
    const map<word, scalar>& speciesCon,
    ChemistryData& data
)
{
    //- Tmp M
    scalar M{0};

    //- ENHANCED species
    wordList enhancedspecies;

    const bool& EH = data.ENHANCED(r);


    //- Use enhance values
    if (EH)
    {
        //- Get information
        enhancedspecies = data.enhancedspecies(r);
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
                forAll(enhancedspecies, e)
                {
                    //- If found, then use modified value
                  //  if (enhancedspecies[e] == species[s])
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


*/

