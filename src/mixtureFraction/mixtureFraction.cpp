/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2015 Holzmann-cfd
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License This file is part of Automatic Flamelet Constructor.

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
    Chemistry& chem,
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
    if (prop.input() == "mol")
    {

        //- Calculate mole fraction (initial linear distribution)
        {
            const wordList& species = chem.species();

            //- Set the oxidizer mol fraction
            map<word, scalar> oxidizerMolFraction = prop.oxidizerCompMol();
            //
            //- Set the fuel mol fraction
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
                Info<< "Discrete point Z (mol): " << Z_ << endl;

                forAll(species, i)
                {
                    if (speciesMol_.at(species[i]) > 0)
                    {
                        Info<< species[i] << ": " 
                            << speciesMol_.at(species[i]) << endl;
                    }
                }
            }

            //- Calculate mass fraction Y
            XtoY();
        }
    }

    // Initialize with mass fraction
    else if (prop.input() == "mass")
    {
        //- Calculate mole fraction (initial linear distribution)
        {
            const wordList& species = chem.species();

            //- Set the oxidizer mass fraction
            map<word, scalar> oxidizerMassFraction = prop.oxidizerCompMass();
            //
            //- Set the fuel mass fraction
            map<word, scalar> fuelMassFraction = prop.fuelCompMass();

            //- Set all species to zero
            forAll(species, i)
            {
                if (oxidizerMassFraction[species[i]] > 1e-10)
                {
                    speciesMass_[species[i]] =
                        oxidizerMassFraction[species[i]]
                      * (-1 * Zvalue + 1);
                }
                else if (fuelMassFraction[species[i]] > 1e-10)
                {
                    speciesMass_[species[i]] =
                        fuelMassFraction[species[i]]
                      * Zvalue;
                }
                else
                {
                    speciesMass_[species[i]] = 0;
                }
            }

            if (debug)
            {
                Info<< "Discrete point Z (mass): " << Z_ << endl;

                forAll(species, i)
                {
                    if (speciesMass_.at(species[i]) > 0)
                    {
                        Info<< species[i] << ": " 
                            << speciesMass_.at(species[i]) << endl;
                    }
                }
            }
            //- Calculate mol fraction X
            YtoX();
        }
    }

    //- Initiate all other stuff
    {
        //- Calculate concentration [X]
        XtoC();

        //- Calculate mean density rho
        rhoX();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::MixtureFraction::~MixtureFraction()
{
    if (debug)
    {
//        Info<< "Destruct MixtureFraction\n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AFC::MixtureFraction::calculateMeanMW
(
    const word& funcSW
)
{
    const wordList& species = chemistry_.species();        

    //- Reset
    MW_ = 0;

    //- Calculate mean moleculare weight with mol fraction
    if (funcSW == "mol")
    {
        forAll(species, s)
        {
            MW_ += speciesMol_.at(species[s]) * thermo_.MW(species[s]);
        }

        updatedMW_ = true; 
    } 
    
    //- Calculate mean moleculare weight with mass fraction
    else if (funcSW == "mass")
    {
        scalar tmp{0};

        forAll(species, s)
        {
            tmp += speciesMass_.at(species[s]) / thermo_.MW(species[s]);
        }

        MW_ = 1/tmp;
    }
    //- Calculate mean moleculare weight with concentration 
    else if (funcSW == "con")
    {
        scalar numerator{0};
        scalar denominator{0};

        forAll(species, s)
        {
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


void AFC::MixtureFraction::calculateMeanCp
(
    const scalar& T 
)
{
    const wordList& species = chemistry_.species();        

    //- Reset
    cp_ = 0;

    forAll(species, s)
    {
        cp_ += speciesMass_.at(species[s]) * thermo_.cp(species[s], T);
    }
}


void AFC::MixtureFraction::calculateMeanH
(
    const scalar& T 
)
{
    const wordList& species = chemistry_.species();        

    //- Reset
    H_ = 0;

    forAll(species, s)
    {
        H_ += speciesMass_.at(species[s]) * thermo_.H(species[s], T);
    }
}


void AFC::MixtureFraction::calculateMeanS
(
    const scalar& T 
)
{
    const wordList& species = chemistry_.species();        

    //- Reset
    S_ = 0;

    forAll(species, s)
    {
        S_ += speciesMass_.at(species[s]) * thermo_.S(species[s], T);
    }
}


void AFC::MixtureFraction::calculateMeanG
(
    const scalar& T 
)
{
    const wordList& species = chemistry_.species();        

    //- Reset
    G_ = 0;

    forAll(species, s)
    {
        G_ += speciesMass_.at(species[s]) * thermo_.G(species[s], T);
    }
}


void AFC::MixtureFraction::calculateMeanG
(
    const scalar& H,
    const scalar& S,
    const scalar& T 
)
{
    //- Calculate Gibbs energy with known values
    G_ = thermo_.G(H, S, T);
}


void AFC::MixtureFraction::updatekfkb
(
    const scalar& T    
)
{
    //- For the update procedure we need thermo and chemistry stuff
    chemistry_.updatekfkb(T, thermo_);
}


void AFC::MixtureFraction::calculateHf
(
    const word& species,
    const scalar& T
) const
{
    Info<< "dHf("<<species<<") = " <<  thermo_.Hf(species, T) * AFC::Constants::jouleToCal / 1000 << endl;
}


// * * * * * * * * * * * * Conversation Functions  * * * * * * * * * * * * * //

void AFC::MixtureFraction::YtoX()
{
    const wordList& species = chemistry_.species();

    //- update mean molecular weight
    calculateMeanMW("mass");

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

        Info<< "    YtoX(), sum of mol = " << sum << endl;
    }
}


void AFC::MixtureFraction::XtoY()
{
    const wordList& species = chemistry_.species();

    // update mean molecular weight
    calculateMeanMW("mol");

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
            if (speciesMass_.at(species[s]) > 0)
            {
                Info<< species[s] << ": " << speciesMass_.at(species[s])
                    << endl;
            }

            sum += speciesMass_.at(species[s]);
        }

        Info<< "    XtoY(), sum of mass = " << sum << endl;
    }
}


void AFC::MixtureFraction::YtoC()
{
    const wordList& species = chemistry_.species();

    // update mean molecular weight
    calculateMeanMW("mass");

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

        Info<< "    YtoCon(), sum of concentration = " << sum << endl;
    }
}


void AFC::MixtureFraction::XtoC()
{
    const wordList& species = chemistry_.species();

    // update mean molecular weight
    calculateMeanMW("mol");

    scalar XT{0};

    forAll(species, s)
    {
        XT += speciesMol_.at(species[s]) * T();
    }

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

        Info<< "    XtoCon(), sum of concentration = " << sum << endl;
    }
}


void AFC::MixtureFraction::rhoX()
{
    const wordList& species = chemistry_.species();

    scalar XT{0};

    forAll(species, s)
    {
        XT += speciesMol_.at(species[s]) * T();
    }

    const scalar& p = properties_.p();

    //- Reset rho_
    rho_ = 0;
   
    forAll(species, s)
    {
       rho_ += speciesMol_.at(species[s]) * p 
           / (AFC::Constants::R * XT) * thermo_.MW(species[s]);
    }

    if (debug)
    {
        Info<< "    Mean density (rhoX): " << rho_ << endl;
    }
}


void AFC::MixtureFraction::rhoY()
{
    const wordList& species = chemistry_.species();

    scalar YTByMW{0};

    forAll(species, s)
    {
        YTByMW += speciesMass_.at(species[s]) * T() / thermo_.MW(species[s]);
    }

    //- Reset rho_
    rho_ = 0;

    const scalar& p = properties_.p();

    forAll(species, s)
    {
        rho_ += (p * speciesMass_.at(species[s])) 
             / (AFC::Constants::R * YTByMW); 
    }

    if (debug)
    {
        Info<< "    Mean density (rhoY): " << rho_ << endl;
    }
}


void AFC::MixtureFraction::rhoC()
{
    const wordList& species = chemistry_.species();

    //- Reset rho_
    rho_ = 0;

    forAll(species, s)
    {
        rho_ += speciesCon_.at(species[s]) * thermo_.MW(species[s]);
    }

    if (debug)
    {
        Info<< "    Mean density (rhoC): " << rho_ << endl;
    }
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

AFC::map<AFC::word, AFC::scalar> AFC::MixtureFraction::mol() const
{
    return speciesMol_;
}


AFC::scalar AFC::MixtureFraction::T() const
{
    return temperature_;
}


AFC::scalar AFC::MixtureFraction::cp() const
{
    return cp_;
}


AFC::scalar AFC::MixtureFraction::calculateCp
(
    const word& species,
    const scalar& T
) const
{
    return thermo_.cp(species, T);
}


AFC::scalar AFC::MixtureFraction::H() const
{
    return H_;
}


AFC::scalar AFC::MixtureFraction::calculateH
(
    const word& species,
    const scalar& T
) const
{
    return thermo_.H(species, T);
}


AFC::scalar AFC::MixtureFraction::S() const
{
    return S_; 
}


AFC::scalar AFC::MixtureFraction::calculateS
(
    const word& species,
    const scalar& T
) const
{
    return thermo_.S(species, T);
}


AFC::scalar AFC::MixtureFraction::G() const
{
    return G_;
}


// ************************************************************************* //
