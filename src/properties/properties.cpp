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

#include "properties.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Properties::Properties
(
    const string& fileName
)
:
    mfPoints_(0),

    vmfPoints_(0),

    TOxidizer_(0),

    TFuel_(0)
{
    PropertiesReader mixFracReader(fileName);

    mixFracReader.read(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Properties::~Properties()
{
    if (debug)
    {
        Info<< "Destruct Properties\n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AFC::Properties::insertMFPoints
(
    const int& mfPoints
)
{
    mfPoints_ = mfPoints;
}


void AFC::Properties::insertVMFPoints
(
    const int& vmfPoints
)
{
    vmfPoints_ = vmfPoints;
}


void AFC::Properties::insertScalarDissipationRates
(
    const scalar& sDR
)
{
    sDR_.push_back(sDR);
}


void AFC::Properties::insertTemperatureOxidizer
(
    const scalar& TOxidizer
)
{
    if (debug)
    {
        Info<< "Oxidizer temperature: " << TOxidizer << endl;
    }

    TOxidizer_ = TOxidizer; 
}


void AFC::Properties::insertTemperatureFuel
(
    const scalar& TFuel
)
{
    if (debug)
    {
        Info<< "Fuel temperature: " << TFuel << endl;
    }

    TFuel_ = TFuel;
}


void AFC::Properties::insertCompositionOxidizerMol
(
    const word& species,
    const scalar& molFraction
)
{
    if (debug)
    {
        Info<< "Oxidizer species: " << species << "  " << molFraction << endl;
    }

    speciesOxidizer_.push_back(species);

    oxidizerCompMol_[species] = molFraction;
}


void AFC::Properties::insertCompositionFuelMol
(
    const word& species,
    const scalar& molFraction
)
{
    if (debug)
    {
        Info<< "Fuel species: " << species << "  " << molFraction << endl;
    }

    speciesFuel_.push_back(species);

    fuelCompMol_[species] = molFraction;
}


// * * * * * * * * * * * * * * * * Check function  * * * * * * * * * * * * * //

void AFC::Properties::check()
{
    if (mfPoints_ == 0)
    {
        FatalError
        (
            "    No mixtureFractionPoints defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if (vmfPoints_ == 0)
    {
        FatalError
        (
            "    No varianzOfMixtureFractionPoints defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if (TOxidizer_ == 0)
    {
        FatalError
        (
            "    No temperature for oxidizer defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if (TFuel_ == 0)
    {
        FatalError
        (
            "    No temperature for fuel defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if (sDR_.empty())
    {
        FatalError
        (
            "    No scalar dissipation rates defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if (oxidizerCompMol_.empty())
    {
        FatalError
        (
            "    No oxidizer species defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if (fuelCompMol_.empty())
    {
        FatalError
        (
            "    No fuel species defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    scalar sum{0};

    forAll(speciesOxidizer_, i)
    {
        sum += oxidizerCompMol_[speciesOxidizer_[i]];
    } 

    scalar epsilon{1e-12};

    //- Double comparison; with epsilon
    if ((1.0 - sum) > epsilon)
    {
        FatalError
        (
            "    Sum of mole fraction of oxidizer is not 1 ("
            + std::to_string(sum) + ").\n"
            "    Please check your fraction composition for the oxidizer "
            "stream in your afcDict.",
            __FILE__,
            __LINE__
        );
    }

    sum = 0;
    forAll(speciesFuel_, i)
    {
        sum += fuelCompMol_[speciesFuel_[i]];
    }

    //- Double comparison; with epsilon
    if ((1.0 - sum) > epsilon)
    {
        FatalError
        (
            "    Sum of mole fraction of fuel is not 1 ("
            + std::to_string(sum) + ").\n"
            "    Please check your fraction composition for the fuel "
            "stream in your afcDict.",
            __FILE__,
            __LINE__
        );
    }
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

AFC::wordList AFC::Properties::speciesOxidizer() const
{
    return speciesOxidizer_;
}


AFC::wordList AFC::Properties::speciesFuel() const
{
    return speciesFuel_;
}


// ************************************************************************* //
