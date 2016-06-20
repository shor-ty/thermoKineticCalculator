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
#include <fstream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::MixtureFraction::MixtureFraction
(
    Chemistry& chem,
    const Thermo& therm,
    const Transport& trans,
    const Properties& prop,
    const scalar& sDR,
    const scalar& defect
)
: 
    defect_(defect),

    sDR_(sDR),

    thermo_(therm),

    transport_(trans),

    chemistry_(chem),

    properties_(prop)
{
    //- nPoints
    const int& nZ = prop.nZPoints();
    const scalar delta = 1. / nZ;

    //- Resize fields
    {
        Z_.resize(nZ+1);
        T_.resize(nZ+1);
        rho_.resize(nZ+1);
        cp_.resize(nZ+1);
        MW_.resize(nZ+1);
        mu_.resize(nZ+1);
        lambda_.resize(nZ+1);
        Cmix_.resize(nZ+1);
    }

    //- Init fields with zeros
    forEach(Z_, i)
    {
        Z_[i] = 0;
        T_[i] = 0;
        rho_[i] = 0;
        cp_[i] = 0;
        MW_[i] = 0;
        mu_[i] = 0;
        lambda_[i] = 0;
        Cmix_[i] = 0;
    }

    //- Resize map fields
    {
        //- Species in chemistry 
        const wordList& species = chem.species();

        //- Generate tmp map with all species and zero value
        map<word, scalar> tmp;

        forAll(species, s)
        {
            tmp[s] = 0; 
        }

        //- Resize and init
        forEach(Z_, i)
        {
            speciesMol_.push_back(tmp);
            speciesMass_.push_back(tmp);
            speciesCon_.push_back(tmp);
        }
    }
        

    //- Calculate values for Z
    for (int i=0; i<=nZ; i++)
    {
        Z_[i] = i * delta;
    }


    //- Calculate temperature profile (linear)
    //  + Calculate concentration of each point
    forEach(T_, i)
    {
        scalar oxidizerTemperature = prop.oxidizerTemperature();
        scalar fuelTemperature = prop.fuelTemperature(); 

        T_[i] =
            (fuelTemperature - oxidizerTemperature)
          * Z_[i] + oxidizerTemperature;
        
        Cmix_[i] = C(T_[i]);
    }


    // Initialize with mole fraction
    /*if (prop.input() == "mol")
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
            forAll(species, s)
            {
                if
                (
                    oxidizerMolFraction[s] > 1e-10
                 || fuelMolFraction[s] > 1e-10
                )
                {
                    speciesMol_[s] =
                        oxidizerMolFraction[s] * (1 - Zvalue)
                      + fuelMolFraction[s] * Zvalue;
                }
                else
                {
                    speciesMol_[s] = 0;
                }
            }

            if (debug)
            {
                Info<< "Discrete point Z (mol): " << Z_ << endl;

                forAll(species, s)
                {
                    if (speciesMol_.at(s) > 0)
                    {
                        Info<< s << ": " << speciesMol_.at(s) << endl;
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
            forAll(species, s)
            {
                if (oxidizerMassFraction[s] > 1e-10)
                {
                    speciesMass_[s] =
                        oxidizerMassFraction[s]
                      * (-1 * Zvalue + 1);
                }
                else if (fuelMassFraction[s] > 1e-10)
                {
                    speciesMass_[s] =
                        fuelMassFraction[s]
                      * Zvalue;
                }
                else
                {
                    speciesMass_[s] = 0;
                }
            }

            if (debug)
            {
                Info<< "Discrete point Z (mass): " << Z_ << endl;

                forAll(species, s)
                {
                    if (speciesMass_.at(s) > 0)
                    {
                        Info<< s << ": " << speciesMass_.at(s) << endl;
                    }
                }
            }
            //- Calculate mol fraction X
            YtoX();
        }
    }*/

    //- Initiate all other stuff
    {
        //- Calculate concentration [X]
//        XtoC();

        //- Calculate mean density rho
//        rhoX();
    }

    //- Generate a log file
    summary();
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

int AFC::MixtureFraction::nZPoints() const
{
    return properties_.nZPoints();
}


AFC::scalar AFC::MixtureFraction::Zvalue
(
    const int& x
)
{
    return Z_[x];
}


AFC::scalarField& AFC::MixtureFraction::T() 
{
    return T_;
}


/*void AFC::MixtureFraction::calculateMeanMW
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


void AFC::MixtureFraction::calculateOmega
(
    const word& species,
    const scalar& T,
    map<word, scalar>& con
)
{
    chemistry_.calculateOmega(T, con, thermo_);
}
*/

/*AFC::scalar AFC::MixtureFraction::calculateHf
(
    const word& species,
    const scalar& T
) const
{
    return thermo_.Hf(species, T);
}*/


// * * * * * * * * * * * * * * * Update Functions  * * * * * * * * * * * * * //

void AFC::MixtureFraction::updateRho()
{
    //- Calculate rho using concentration [X]
    rhoC();
}


void AFC::MixtureFraction::updateCp()
{
    //- Species
    const wordList species = chemistry_.species();

    //- Reset cp
    Cp_ = 0;

    forAll(species, s)
    {
        Cp_ += calculateCp(species[s], temperature_)
            * speciesMol_.at(species[s]); 
    }
}

void AFC::MixtureFraction::updateC()
{
    YtoC();
}


// * * * * * * * * * * * * Conversation Functions  * * * * * * * * * * * * * //

/*void AFC::MixtureFraction::YtoX()
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

    //- Pressure [Pa]
    const scalar& p = properties_.p();

    //- Unit of concentration [mol/m^3]
    //  Change unit to [mol/cm^3] factor 100cm * 100cm * 100cm / m^3
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
}*/


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

AFC::scalar AFC::MixtureFraction::C
(
    const scalar& T 
) const
{
    return thermo_.C(T);
}


/*AFC::map<AFC::word, AFC::scalar> AFC::MixtureFraction::mol() const
{
    return speciesMol_;
}


AFC::map<AFC::word, AFC::scalar> AFC::MixtureFraction::mol() const
{
    return speciesMol_;
}*/


AFC::map<AFC::word, AFC::scalar>& AFC::MixtureFraction::mass()
{
    return speciesMass_;
}


AFC::map<AFC::word, AFC::scalar> AFC::MixtureFraction::con()
{
    return speciesCon_;
}


AFC::scalar& AFC::MixtureFraction::T()
{
    return temperature_;
}


AFC::scalar AFC::MixtureFraction::T() const
{
    return temperature_;
}


AFC::scalar AFC::MixtureFraction::cp() const
{
    return cp_;
}


AFC::scalar AFC::MixtureFraction::Cp() const
{
    return Cp_;
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
*/


// * * * * * * * * * * * * * * Summary function  * * * * * * * * * * * * * * //

void AFC::MixtureFraction::summary() const
{
    std::filebuf file;

    //- Open file
    file.open("analyze", std::ios::out);

    std::ostream data(&file);

    //- Set scientific notation
    data.setf(std::ios::scientific, std::ios::floatfield);

    //- Header
    data<< Header() << "\n"; 

    //- Build the chemistry summary
    chemistry_.summary(data);

    //- Build the thermodynamic summary
    thermo_.summary(data);


    //file<< dataOfChemistry;

    file.close();
}

// ************************************************************************* //
