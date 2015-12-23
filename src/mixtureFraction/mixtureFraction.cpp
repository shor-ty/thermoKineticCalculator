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
#include "constants.hpp"

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

    Z_(Zvalue),

    thermo_(therm),

    transport_(trans),

    chemistry_(chem),

    properties_(prop)
{
    //- Calculate temperature profile (linear)
    {
        scalar oxidizerTemperature = prop.oxidizerTemperature();
        scalar fuelTemperature = prop.fuelTemperature(); 

        temperature_ =
            (fuelTemperature - oxidizerTemperature)
          * Zvalue
          + oxidizerTemperature;
    }

    // Initialize with mole fraction
    // if ()
    // {

        //- Calculate mole fraction (initial linear distribution)
        {
            const wordList& species = chem.species();

            //- Set the oxidizer mole fraction
            map<word, scalar> oxidizerMolFraction = prop.oxidizerCompMol();
            //
            //- Set the fuel mole fraction
            map<word, scalar> fuelMolFraction = prop.fuelCompMol();

            //- Set all species to zero
            forAll(species, i)
            {
                if (oxidizerMolFraction[species[i]] > 1e-10)
                {
                    speciesMol_[species[i]] =
                        oxidizerMolFraction[species[i]]
                      * (-1 * Zvalue + 1);
                }
                else if (fuelMolFraction[species[i]] > 1e-10)
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

            if (debug)
            {
                Info<< "Discrete point Z: " << Z_ << endl;

                forAll(species, i)
                {
                    if (speciesMol_.at(species[i]) > 0)
                    {
                        Info<< species[i] << ": " << speciesMol_.at(species[i]) << endl;
                    }
                }
            }

            //- Calculate mass fraction Y
            XtoY();

            //- Calculate concentration [X]
            XtoC();

            //- Calculate mean density rho
            XtoRho();

            if (debug)
            {
                Info<< endl;
            }
        }
    //}

    // Initialize with mass fraction
    // else if ()
    // {
    // }
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

void AFC::MixtureFraction::calcMeanMW
(
    const word& funcSW
)
{
    const wordList& species = chemistry_.species();        

    //- Calculate mean moleculare weight with mol fraction
    if (funcSW == "mol")
    {
        forAll(species, s)
        {
            if (debugMW)
            {
                Info<< "MW calc, Species: " << species[s]
                    << " << MW: " << thermo_.MW(species[s])
                    << " << X: " << speciesMol_.at(species[s]) << endl;
            }

            MW_ += speciesMol_.at(species[s])*thermo_.MW(species[s]);
        }
    }
    //- Calculate mean moleculare weight with mass fraction
    else if (funcSW == "mass")
    {
        scalar tmp{0};

        forAll(species, s)
        {
            if (debugMW)
            {
                Info<< "MW calc, Species: " << species[s]
                    << " << MW: " << thermo_.MW(species[s])
                    << " << Y: " << speciesMass_.at(species[s]) << endl;
            }

            tmp += speciesMass_.at(species[s])*thermo_.MW(species[s]);
        }

        MW_ = 1/tmp;
    }
    else if (funcSW == "con")
    {
        scalar numerator{0};
        scalar denominator{0};

        forAll(species, s)
        {
            if (debugMW)
            {
                Info<< "MW calc, Species: " << species[s]
                    << " << MW: " << thermo_.MW(species[s])
                    << " << [X]: " << speciesCon_.at(species[s]) << endl;
            }

            numerator += speciesCon_.at(species[s])*thermo_.MW(species[s]);
            denominator += speciesCon_.at(species[s]);
        }

        MW_ = numerator / denominator;
    }
    else
    {
        FatalError
        (
            "    Calculation of mean moleculare weight not implemented.\n"
            "    You can either use 'mass', 'mol' or 'con' as argument.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::MixtureFraction::YtoX()
{
    const wordList& species = chemistry_.species();

    //- update mean molecular weight
    calcMeanMW("mass");

    forAll(species, s)
    {
        speciesMol_[species[s]]
            = speciesMass_.at(species[s]) * MW_ / thermo_.MW(species[s]);
    }

    if (debug)
    {
        scalar sum{0};

        forAll(species, s)
        {
            sum += speciesMol_.at(species[s]);
        }

        Info<< "YtoX(), sum of mol = " << sum << endl;
    }
}


void AFC::MixtureFraction::XtoY()
{
    const wordList& species = chemistry_.species();

    // update mean molecular weight
    calcMeanMW("mol");

    forAll(species, s)
    {
        speciesMass_[species[s]]
            = speciesMol_.at(species[s]) * thermo_.MW(species[s]) / MW_;
    }

    if (debug)
    {
        scalar sum{0};

        forAll(species, s)
        {
            sum += speciesMass_.at(species[s]);
        }

        Info<< "XtoY(), sum of mass = " << sum << endl;
    }
}


void AFC::MixtureFraction::YtoC()
{
    const wordList& species = chemistry_.species();

    // update mean molecular weight
    calcMeanMW("mass");

    scalar YTMW{0};

    forAll(species, s)
    {
        YTMW += speciesMass_.at(species[s]) * T() / thermo_.MW(species[s]);
    }

    const scalar& p = properties_.p();

    forAll(species, s)
    {
        speciesCon_[species[s]]
            = (p * speciesMass_.at(species[s]) / thermo_.MW(species[s]))
            / (AFC::Constants::R * YTMW);
    }

    if (debug)
    {
        scalar sum{0};

        forAll(species, s)
        {
            sum += speciesCon_.at(species[s]);
        }

        Info<< "YtoCon(), sum of concentration = " << sum << endl;
    }
}


void AFC::MixtureFraction::XtoC()
{
    const wordList& species = chemistry_.species();

    // update mean molecular weight
    calcMeanMW("mol");

    scalar XT{0};

    forAll(species, s)
    {
        XT += speciesMol_.at(species[s]) * T();
    }

    Info<< "XT: " << XT << endl;

    const scalar& p = properties_.p();

    forAll(species, s)
    {
        speciesCon_[species[s]] 
            = speciesMol_.at(species[s]) * p / (XT * AFC::Constants::R);
    }

    if (debug)
    {
        scalar sum{0};

        forAll(species, s)
        {
            sum += speciesCon_.at(species[s]);
        }

        Info<< "XtoCon(), sum of concentration = " << sum << endl;
    }
}


void AFC::MixtureFraction::XtoRho()
{
    const wordList& species = chemistry_.species();

    scalar XT{0};

    forAll(species, s)
    {
        XT += speciesMol_.at(species[s]) * T();
    }

    scalar tmp{0};

    const scalar& p = properties_.p();
   
    forAll(species, s)
    {
       tmp += speciesMol_.at(species[s]) * p 
           / (AFC::Constants::R * XT) * thermo_.MW(species[s]);
    }

    rho_ = tmp;

    if (debug)
    {
        Info<< "   Mean density: " << rho_ << endl;
    }
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

void AFC::MixtureFraction::mols(const word spec) const
{
//    Info<< "Species " << spec << ": " << speciesMol_.at(spec) << endl;   
}


/*AFC::scalar AFC::MixtureFraction::mol
(
    const word& species
)
{
    return speciesMol_[species];
}

*/
AFC::map<AFC::word, AFC::scalar> AFC::MixtureFraction::mol() const
{
    return speciesMol_;
}


AFC::scalar AFC::MixtureFraction::T() const
{
    return temperature_;
}


// ************************************************************************* //
