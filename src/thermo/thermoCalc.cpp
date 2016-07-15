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

#include "thermoCalc.hpp"
#include "constants.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::ThermoCalc::ThermoCalc()
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::ThermoCalc::~ThermoCalc()
{
}


// * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * * //

AFC::scalar AFC::ThermoCalc::MmeanX
(
    const map<word, scalar>& X,
    const map<word, scalar>& MW
) const
{
    scalar tmp{0};

    forAll(X, s)
    {
        tmp += s.second * MW.at(s.first); 
    }

    return tmp;
}


AFC::scalar AFC::ThermoCalc::MmeanY
(
    const map<word, scalar>& Y,
    const map<word, scalar>& MW
) const
{
    scalar tmp{0};

    forAll(Y, s)
    {
        tmp += s.second / MW.at(s.first); 
    }

    return pow(tmp, -1);
}


AFC::scalar AFC::ThermoCalc::MmeanC
(
    const map<word, scalar>& C,
    const map<word, scalar>& MW
) const
{
    scalar tmp{0};
    scalar sumC{0};

    forAll(C, s)
    {
        tmp += s.second * MW.at(s.first); 
        sumC += s.second;
    }

    return tmp/sumC;
}


AFC::scalar AFC::ThermoCalc::rhoMean
(
    const scalar& Mmean,
    const scalar& p,
    const scalar& T
) const
{
    return (p * Mmean/ (AFC::Constants::R * T));
}


AFC::scalar AFC::ThermoCalc::C
(
    const scalar& p,
    const scalar& T
) const
{
    return (p / (AFC::Constants::R * T));
}


AFC::scalar AFC::ThermoCalc::cp
(
    const word& species,
    const scalar& T,
    const ThermoData& data
) const
{
    scalarField coeffs = getCoeffs(species, T, data);

    //- calculate and return [J/mol/K]
    return
    (
        (
            coeffs[0]
          + coeffs[1] * T
          + coeffs[2] * pow(T, 2)
          + coeffs[3] * pow(T, 3)
          + coeffs[4] * pow(T, 4)
        ) * AFC::Constants::R
    );
}


AFC::scalar AFC::ThermoCalc::H
(
    const word& species,
    const scalar& T,
    const ThermoData& thermoData
) const
{
    scalarField coeffs = getCoeffs(species, T, thermoData);

    //- calculate and return [J/mol]
    return
    (
        (
            coeffs[0]
          + coeffs[1] * T / 2
          + coeffs[2] * pow(T, 2) / 3
          + coeffs[3] * pow(T, 3) / 4
          + coeffs[4] * pow(T, 4) / 5
          + coeffs[5] / T
        ) * AFC::Constants::R * T
    );
}


AFC::scalar AFC::ThermoCalc::dH
(
    const word& species,
    const scalar& T,
    const ThermoData& data 
) const
{
    return (H(species, T, data) - Hf(species, data));
}


AFC::scalar AFC::ThermoCalc::S
(
    const word& species,
    const scalar& T,
    const ThermoData& thermoData
) const
{
    const scalarField coeffs = getCoeffs(species, T, thermoData);

    //- calculate and return [J/mol/K]
    return
    (
        (
            coeffs[0] * log(T)
          + coeffs[1] * T
          + coeffs[2] * pow(T, 2) / 2
          + coeffs[3] * pow(T, 3) / 3
          + coeffs[4] * pow(T, 4) / 4
          + coeffs[6]
          + log(thermoData.p()/AFC::Constants::p0)
        ) * AFC::Constants::R
    );
}


AFC::scalar AFC::ThermoCalc::G
(
    const word& species,
    const scalar& T,
    const ThermoData& thermoData
) const
{
    //- Enthalpy
    const scalar h = H(species, T, thermoData);

    //- Entropy
    const scalar s = S(species, T, thermoData);

    //- calculate free Gibbs energy
    return ( h - s * T);
}


AFC::scalar AFC::ThermoCalc::dG
(
    const word& species,
    const scalar& T,
    const ThermoData& data
) const
{
    return (G(species, T, data) - Gf(species, data));
}


AFC::scalar AFC::ThermoCalc::G
(
    const scalar& h,
    const scalar& s,
    const scalar& T
) const
{
    //- calculate mean free Gibbs energy
    return ( h - s * T);
}


AFC::scalar AFC::ThermoCalc::Hf
(
    const word& species,
    const ThermoData& data
) const
{
    //- Get information about species (atoms)
/*    const wordList& atoms = data.speciesAtoms(species);

    //- Get information about multiplicator
    const scalarList& factors = data.atomFactors(species);

    //- 
    const map<word, scalar>& factor = data.atomFactors();

    //- New multiplicator based on stable species
    scalarList newFactors;

    //- New species word list
    wordList newAtoms;

    //- First: Stable atoms O2, H2, N2 etc...
    forAll(atoms, a)
    {
        //- FACTORS:: If we find these atoms we divide factors by 2
        {
            if
            (
                a == "O"
             || a == "H"
             || a == "N"
            )
            {
                newFactors.push_back(factor.at(a)/2);
            }
            else if (a == "AR")
            {
                newFactors.push_back(factor.at(a));
            }
            else
            {
                newFactors.push_back(factor.at(a));
            }
        }
        //- ATOMS:: Convert to stable atoms
        {
            if (a == "O")
            {
                newAtoms.push_back("O2");
            }
            else if (a == "H")
            {
                newAtoms.push_back("H2");
            }
            else if (a == "N")
            {
                newAtoms.push_back("N2");
            }
            else if (a == "C")
            {
                newAtoms.push_back("CSOLID");
            }
            else
            {
                newAtoms.push_back(a);
            }
        }
    }

    forAll(newFactors, a)
    {
//        Info<< atoms[a] << " -> " << factors[a] << " | new -> " << newFactors[a] << endl;
    }

    //- Second: calculate temperature depended H of species
    const scalar Hspecies = H(species, T, data);

    //- Third: calculate temperature dependend H of stable atoms 
    //  + H2
    //  + O2
    //  + N2
    scalarList Hatoms;

    forAll(newAtoms, a)
    {
        Hatoms.push_back(H(a, T, data));
    }

    scalar deltaHf = Hspecies;

    //- Fourth: Calculate delta Hf of species
    forAll(newAtoms, a)
    {
       deltaHf -= Hatoms[a] * newFactors[a]; 
    }

    return deltaHf;*/
    return H(species, 298, data);
}


AFC::scalar AFC::ThermoCalc::Gf
(
    const word& species,
    const ThermoData& data
) const
{
    return G(species, 298, data);
}


// * * * * * * * * * * * * * * *  Return Functions * * * * * * * * * * * * * //

AFC::scalarField AFC::ThermoCalc::getCoeffs
(
    const word& species,
    const scalar& T,
    const ThermoData& thermoData
) const
{
    //- get temperature range to choose NASA polynomials   
    //  + true -> high temp
    //  + false -> low temp
    const bool highTemp = whichTempRange(species, T, thermoData);

    if (highTemp)
    {
        return thermoData.NASACoeffsHT(species);
    }
    else
    {
        return thermoData.NASACoeffsLT(species);
    }
}


bool AFC::ThermoCalc::whichTempRange
(
    const word& species,
    const scalar& T,
    const ThermoData& thermoData
) const
{
    //- Low temperature 
    const scalar& TL = thermoData.LT(species);

    //- Common temperature
    const scalar& TC = thermoData.CT(species);

    //- High temperature
    const scalar& TH = thermoData.HT(species);

    //- Check if temperature T is in range of NASA
    if (T < TL)
    {
        if (warnings)
        {
            Warning
            (
                "   Temperature is lower than the range of NASA polynomials.",
                __FILE__,
                __LINE__
            );
        }

        //- Lower than TL
        return false;
    }
    else if (T > TH)
    {
        if (warnings)
        {
            Warning
            (
                "   Temperature is higher than the range of NASA polynomials.",
                __FILE__,
                __LINE__
            );
        }

        //- Higher than TH
        return true;
    }
    else
    {
        if
        (
            T >= TL
         && T <= TC
        )
        {
            //- TL to TC 
            return false;
        }
        else
        {
            //- TC to TH
            return true;
        }
    }
}

// ************************************************************************* //
