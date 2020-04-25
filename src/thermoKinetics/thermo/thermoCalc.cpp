/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     T hermo       | Open Source Thermo-Kinetic Library
  c-o-o-c     K iknetic     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
  c     c-o-o-o             |
-------------------------------------------------------------------------------
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

#include "thermoCalc.hpp"
#include "constants.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::ThermoCalc::ThermoCalc
(
    const string fileName,
    const bool thermoInChemistry
)
:
    ThermoData(fileName)
{
    if (debug_)
    {
        Info<< "ThermoCalc Constructor\n" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::ThermoCalc::~ThermoCalc()
{
    if (debug_)
    {
        Info<< "ThermoCalc Destructor\n" << endl;
    }
}


// * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * * //

TKC::scalar TKC::ThermoCalc::MWmeanX
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


TKC::scalar TKC::ThermoCalc::MWmeanY
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


TKC::scalar TKC::ThermoCalc::MWmeanC
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


TKC::scalar
TKC::ThermoCalc::rho(const scalar T, const scalar p, const scalar MW) const
{
    return (p * MW / TKC::Constants::R / T);
}


TKC::scalar TKC::ThermoCalc::rhoMean
(
    const scalar Mmean,
    const scalar p,
    const scalar T
) const
{
    return (p * Mmean/ (TKC::Constants::R * T));
}


TKC::scalar TKC::ThermoCalc::C(const scalar p, const scalar T) const
{
    return (p / (TKC::Constants::R * T));
}


TKC::scalar TKC::ThermoCalc::cp(const word species, const scalar T) const
{
    const scalarField& coeffs = getCoeffs(species, T);

    //- calculate and return [J/mol/K]
    return
    (
        (
            coeffs[0]
          + coeffs[1] * T
          + coeffs[2] * pow(T, 2)
          + coeffs[3] * pow(T, 3)
          + coeffs[4] * pow(T, 4)
        ) * TKC::Constants::R
    );
}


TKC::scalar TKC::ThermoCalc::cv(const word species, const scalar T) const
{
    const scalarField& coeffs = getCoeffs(species, T);

    //- calculate and return [J/mol/K]
    return
    (
        (
            coeffs[0]
          + coeffs[1] * T
          + coeffs[2] * pow(T, 2)
          + coeffs[3] * pow(T, 3)
          + coeffs[4] * pow(T, 4)
        ) * TKC::Constants::R
        - TKC::Constants::R
    );
}


TKC::scalar TKC::ThermoCalc::h(const word species, const scalar T) const
{
    const scalarField& coeffs = getCoeffs(species, T);

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
        ) * TKC::Constants::R * T
    );
}


TKC::scalar TKC::ThermoCalc::dhf(const word species, const scalar T) const
{
    return (h(species, T) - hf(species));
}


TKC::scalar TKC::ThermoCalc::s(const word species, const scalar T) const
{
    const scalarField& coeffs = getCoeffs(species, T);

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
          + log(p()/TKC::Constants::p0)
        ) * TKC::Constants::R
    );
}


TKC::scalar TKC::ThermoCalc::g(const word species, const scalar T) const
{
    //- Enthalpy
    const scalar th = h(species, T);

    //- Entropy
    const scalar ts = s(species, T);

    //- calculate free Gibbs energy
    return ( th - ts * T);
}


TKC::scalar TKC::ThermoCalc::dgf(const word species, const scalar T) const
{
    return (g(species, T) - gf(species));
}


TKC::scalar
TKC::ThermoCalc::g(const scalar h, const scalar s, const scalar T) const
{
    //- calculate mean free Gibbs energy
    return ( h - s * T);
}


TKC::scalar TKC::ThermoCalc::hf(const word species) const
{
    return h(species, 298);
}


TKC::scalar TKC::ThermoCalc::gf(const word species) const
{
    return g(species, 298);
}


TKC::scalar TKC::ThermoCalc::h0(const word species, const scalar T) const
{
    const scalarField& coeffs = getCoeffs(species, T);

    //- calculate and return [J/mol/K]
    return
    (
        coeffs[6] * TKC::Constants::R
    );
}


// * * * * * * * * * * * * * * *  Return Functions * * * * * * * * * * * * * //

TKC::scalarField TKC::ThermoCalc::getCoeffs
(
    const word species,
    const scalar T
) const
{
    //- get temperature range to choose NASA polynomials   
    //  + true -> high temp
    //  + false -> low temp
    const bool highTemp = whichTempRange(species, T);

    if (highTemp)
    {
        return NASACoeffsHT(species);
    }
    else
    {
        return NASACoeffsLT(species);
    }
}


bool TKC::ThermoCalc::whichTempRange(const word species, const scalar T) const
{
    //- Low temperature 
    const scalar TL = LT(species);

    //- Common temperature
    const scalar TC = CT(species);

    //- High temperature
    const scalar TH = HT(species);

    //- Check if temperature T is in range of NASA
    if (T < TL)
    {
        if (warnings_)
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
        if (warnings_)
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
        if (T >= TL && T <= TC)
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
