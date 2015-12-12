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

#include "mixtureFraction.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::MixtureFraction::MixtureFraction
(
    const Chemistry& chem,
    const Thermo& therm,
    const Transport& trans,
    const Properties& prop,
    const scalar& Zvalue,
    const scalar& defect
)
: 
    defect_(defect),

    thermo_(therm),

    transport_(trans)
{

    //- Calculate mole fraction (initial linear distribution)
    {
        wordList species = chem.species();

        //- Set all species to zero
        forAll(species, i)
        {
            speciesMol_[species[i]] = 0;
        }

        //- Set the oxidizer mole fraction
        map<word, scalar> oxidizerMolFraction = prop.oxidizerCompMol();

        forAll(species, i)
        {
            if (oxidizerMolFraction[species[i]] > 1e-10)
            {
                speciesMol_[species[i]] =
                    oxidizerMolFraction[species[i]]
                  * (-1 * Zvalue + 1);
            }
            else
            {
                speciesMol_[species[i]] = 0;
            }
        }

        //- Set the fuel mole fraction
        map<word, scalar> fuelMolFraction = prop.fuelCompMol();

        forAll(species, i)
        {
            if (fuelMolFraction[species[i]] > 1e-10)
            {
                speciesMol_[species[i]] =
                    fuelMolFraction[species[i]]
                  * Zvalue;
            }
            else
            {
                speciesMol_[species[i]] = 0;
            }
        }
    }

    //- Calculate temperature profile (linear)
    {
        scalar oxidizerTemperature = prop.oxidizerTemperature();
        scalar fuelTemperature = prop.fuelTemperature(); 

        temperature_ =
            (fuelTemperature - oxidizerTemperature)
          * Zvalue
          + oxidizerTemperature;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::MixtureFraction::~MixtureFraction()
{
    if (debug)
    {
        Info<< "Destruct MixtureFraction\n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

AFC::scalar AFC::MixtureFraction::mol
(
    const word& species
)
{
    return speciesMol_[species];
}


AFC::scalar AFC::MixtureFraction::T() const
{
    return temperature_;
}


// ************************************************************************* //
