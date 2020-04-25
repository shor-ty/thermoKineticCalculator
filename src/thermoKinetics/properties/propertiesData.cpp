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

#include "propertiesData.hpp"
#include "propertiesReader.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::PropertiesData::PropertiesData
(
    const string fileName
)
{
    PropertiesReader reader(fileName);

    reader.read(*this);

    //- Convert Y to X or X to Y
    /*convertFractions();

    //- Extract the single atomic elements
    {
        oxidizerA_ =
            PropertiesDataCalc::elementDecomposition(speciesOxidizer_, thermo_);


        fuelA_ = PropertiesDataCalc::elementDecomposition(speciesFuel_, thermo_);
    }

    //- Calc the element mass fraction
    {
        oxidizerZj_ =
            PropertiesDataCalc::elementMassFraction
            (
                oxidizerA_,
                oxidizerY_,
                thermo_
            );

        fuelZj_ =
            PropertiesDataCalc::elementMassFraction
            (
                fuelA_,
                fuelY_,
                thermo_
            );


        //- Assign the element mass fractions of the fuel and oxidizer species
        //  Used for the stochiometric calculation

    }

    //- Calc the element mole fraction
    {
        oxidizerWj_ =
            PropertiesDataCalc::elementMolFraction(oxidizerZj_, thermo_);

        fuelWj_ =
            PropertiesDataCalc::elementMolFraction(fuelZj_, thermo_);
    }

    //- Calc the stochiometric fraction Zst
    scalar oxidizer, fuel;

        oxidizer = oxidizerY_.at(oxidizer_);
        fuel = fuelY_.at(fuel_);

        Info<< oxidizer << ":" << fuel << endl;

    Zst_ = PropertiesDataCalc::Zst
        (
            fuelZj_,
            fuelY_.at(fuel_),
            oxidizerY_.at(oxidizer_),
            thermo_
         );

    Info << "Zst = " << Zst_ << endl;


    std::terminate();

    //- Initialize
    //scalar a = PropertiesDataCalc::Zst(fuelY(), oxidizerY());
    */
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::PropertiesData::~PropertiesData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void TKC::PropertiesData::insertInertSpecies(const word inert)
{
    inert_ = inert;
}


void TKC::PropertiesData::insertC
(
    const word species,
    const scalar concentration
)
{
    n_[species] = concentration;
    X_[species] = 0;
    Y_[species] = 0;
}


void TKC::PropertiesData::insertX
(
    const word species,
    const scalar mole
)
{
    n_[species] = 0;
    X_[species] = mole;
    Y_[species] = 0;
}


void TKC::PropertiesData::insertY
(
    const word species,
    const scalar mass
)
{
    n_[species] = 0;
    X_[species] = 0;
    Y_[species] = mass;
}


void TKC::PropertiesData::inputMode(const word mode)
{
    if (mode == "concentration")
    {
        inputConcentration_ = true;
    }
    else if (mode == "mole")
    {
        inputMole_ = true;
    }
    else if (mode == "mass")
    {
        inputMass_ = true;
    }
}


// * * * * * * * * * * * * * * * Other functions * * * * * * * * * * * * * * //


/*void TKC::PropertiesData::check()
{

    scalar sum{0};

    scalar epsilon{1e-12};

    if (inputMol_)
    {
        //- Oxidizer
        forAll(speciesOxidizer_, s)
        {
            sum += oxidizerX_.at(s);
        }

        //- Double comparison; with epsilon
        if ((1.0 - sum) > epsilon)
        {
            ErrorMsg
            (
                "Sum of mol fraction of oxidizer is not 1 ("
                + std::to_string(sum) + ").\n    "
                "Please check your fraction composition for the oxidizer "
                "stream in your afcDict.",
                __FILE__,
                __LINE__
            );
        }

        //- Fuel
        sum = 0;

        forAll(speciesFuel_, s)
        {
            sum += fuelX_.at(s);
        }

        //- Double comparison; with epsilon
        if ((1.0 - sum) > epsilon)
        {
            ErrorMsg
            (
                "Sum of mol fraction of fuel is not 1 ("
                + std::to_string(sum) + ").\n    "
                "Please check your fraction composition for the fuel "
                "stream in your afcDict.",
                __FILE__,
                __LINE__
            );
        }
    }
    else if (inputMass_)
    {
        //- Oxidizer
        forAll(speciesOxidizer_, s)
        {
            sum += oxidizerY_.at(s);
        }

        //- Double comparison; with epsilon
        if ((1.0 - sum) > epsilon)
        {
            ErrorMsg
            (
                "Sum of mass fraction of oxidizer is not 1 ("
                + std::to_string(sum) + ").\n    "
                "Please check your fraction composition for the oxidizer "
                "stream in your afcDict.",
                __FILE__,
                __LINE__
            );
        }

        //-Fuel
        sum = 0;

        forAll(speciesFuel_, s)
        {
            sum += fuelY_.at(s);
        }

        //- Double comparison; with epsilon
        if ((1.0 - sum) > epsilon)
        {
            ErrorMsg
            (
                "Sum of mass fraction of fuel is not 1 ("
                + std::to_string(sum) + ").\n    "
                "Please check your fraction composition for the fuel "
                "stream in your afcDict.",
                __FILE__,
                __LINE__
            );
        }
    }
    else
    {
        ErrorMsg
        (
            "Either mol fraction or mass fraction defined for oxidizer"
            " stream.",
            __FILE__,
            __LINE__
        );
    }

    //- Check if pressure is set
    if (p_ <= 0)
    {
        ErrorMsg
        (
            "Pressure is not set or set not correct.\n    "
            "Please set the pressure in afcDict.",
            __FILE__,
            __LINE__
        );
    }

    //- Check if oxidizer is set and inside the oxidizer species list
    {
        if (oxidizer_.empty())
        {
            ErrorMsg
            (
                "    No oxidizer species is set.\n"
                "    Please check the oxidizer keyword in afcDict.",
                __FILE__,
                __LINE__
            );
        }

        if(oxidizerY_.find(oxidizer_) == oxidizerY_.end())
        {
            ErrorMsg
            (
                "    The oxidizer species '" + oxidizer_ + "' was not found"
                " in the oxidizer mass- or moleFraction dictionary.\n"
                "    Please check the afcDict.",
                __FILE__,
                __LINE__
            );
        }
    }

    //- Check if fuel is set and inside the fuel species list
    {
        if (fuel_.empty())
        {
            ErrorMsg
            (
                "    No fuel species is set.\n"
                "    Please check the fuel keyword in afcDict.",
                __FILE__,
                __LINE__
            );
        }

        if(fuelY_.find(fuel_) == fuelY_.end())
        {
            ErrorMsg
            (
                "    The fuel species '" + fuel_ + "' was not found"
                " in the fuel mass- or moleFraction dictionary.\n"
                "    Please check the afcDict.",
                __FILE__,
                __LINE__
            );
        }
    }
}
*/


void TKC::PropertiesData::convertFractions()
{
    //- Convert to mass fraction
    if (inputMol_)
    {
        XtoY(oxidizerX_, oxidizerY_);
        XtoY(fuelX_, fuelY_);
    }
    //- Convert to mol fraction
    else if (inputMass_)
    {
        YtoX(oxidizerY_, oxidizerX_);
        YtoX(fuelY_, fuelX_);
    }
    //- Already checked within the check() function. However,
    //  maybe something might go wrong - an additional check is okay here
    else
    {
        ErrorMsg
        (
            "Some problem occurs regarding input mode (mass, mole)",
            __FILE__,
            __LINE__
        );
    }
}

void TKC::PropertiesData::XtoY(const map<word, scalar>& X, map<word, scalar>& Y )
{
    scalar M{0};

    //- Calculate M
    loopMap(species, value, X)
    {
       M += value * thermo_.MW(species);
    }

    //- Calculate Y
    loopMap(species, value, X)
    {
        Y.at(species) = value*thermo_.MW(species) / M;
    }
}


void TKC::PropertiesData::YtoX(const map<word, scalar>& Y, map<word, scalar>& X)
{
    scalar M{0};
    scalar YbyM{0};

    //- Calculate YbyM
    loopMap(species, value, Y)
    {
        YbyM += value / thermo_.MW(species);
    }

    M = pow(YbyM, -1);

    loopMap(species, value, Y)
    {
        X.at(species) = M * value / thermo_.MW(species);
    }
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

const TKC::word TKC::PropertiesData::inertSpecies() const
{
    return inertSpecies_;
}


const TKC::map<TKC::word, TKC::scalar> TKC::PropertiesData::C() const
{
    return C_;
}


const TKC::map<TKC::word, TKC::scalar> TKC::PropertiesData::X() const
{
    return X_;
}


const TKC::map<TKC::word, TKC::scalar> TKC::PropertiesData::Y() const
{
    return Y_;
}


const TKC::scalar TKC::PropertiesData::T() const
{
    return T_;
}


const TKC::scalar TKC::PropertiesData::p() const
{
    return p_;
}


const TKC::word TKC::PropertiesData::inputMode() const
{
    if (!inputMol_ && !inputMass_ && !inputConcentration_)
    {
        ErrorMsg
        (
            "Input mode problems. Nor mass, mole or concentration input...",
            __FILE__,
            __LINE__
        );
    }

    word input{"none"};

    if (inputMol_)
    {
        input = "mol";
    }
    else if (inputMass_)
    {
        input = "mass";
    }
    else if (inputConcentration_)
    {
        input = "concentration";
    }

    return input;
}


// ************************************************************************* //
