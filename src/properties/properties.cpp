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
    const string& fileName,
    Thermo& thermo,
    const Chemistry& chemistry
)
:
    thermo_(thermo),

    chemistry_(chemistry)
{
    if (debug_)
    {
        Info<< "Properties Constructor\n" << endl;
    }

    PropertiesReader mixFracReader(fileName);

    mixFracReader.read(*this);

    //- Add pressure to thermoData for better handling
    thermo.p(this->p());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Properties::~Properties()
{
    if (debug_)
    {
        Info<< "Destruct Properties\n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AFC::Properties::insertMFPoints
(
    const int& nZPoints
)
{
    nZPoints_ = nZPoints;
}


void AFC::Properties::insertVMFPoints
(
    const int& nZvarPoints
)
{
    nZvarPoints_ = nZvarPoints;
}


void AFC::Properties::insertEnthalpyDefects
(
    const scalar& defect
)
{
    defects_.push_back(defect);
}


void AFC::Properties::insertScalarDissipationRates
(
    const scalar& sDR
)
{
    sDRs_.push_back(sDR);
}


void AFC::Properties::insertTemperatureOxidizer
(
    const scalar& TOxidizer
)
{
    if (debug_)
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
    if (debug_)
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
    if (debug_)
    {
        Info<< "Oxidizer species (mol): " << species << "  " << molFraction << endl;
    }

    speciesOxidizer_.push_back(species);

    oxidizerX_[species] = molFraction;
}


void AFC::Properties::insertCompositionOxidizerMass
(
    const word& species,
    const scalar& massFraction
)
{
    if (debug_)
    {
        Info<< "Oxidizer species (mass): " << species << "  " << massFraction << endl;
    }

    speciesOxidizer_.push_back(species);

    oxidizerY_[species] = massFraction;
}


void AFC::Properties::insertCompositionFuelMol
(
    const word& species,
    const scalar& molFraction
)
{
    if (debug_)
    {
        Info<< "Fuel species (mol): " << species << "  " << molFraction << endl;
    }

    speciesFuel_.push_back(species);

    fuelX_[species] = molFraction;
}


void AFC::Properties::insertCompositionFuelMass
(
    const word& species,
    const scalar& massFraction
)
{
    if (debug_)
    {
        Info<< "Fuel species (mass): " << species << "  " << massFraction << endl;
    }

    speciesFuel_.push_back(species);

    fuelY_[species] = massFraction;
}


void AFC::Properties::insertRunTime
(
    const scalar& runTime 
)
{
    if (debug_)
    {
        Info<< "Run time of calculation: " << runTime << endl;
    }

    runTime_ = runTime;
}


void AFC::Properties::insertWriteControl
(
    const word& writeControl
)
{
    if (debug_)
    {
        Info<< "Write control: " << writeControl << endl;
    }

    writeControl_ = writeControl;
}


void AFC::Properties::insertWriteControlInterval
(
    const scalar& writeControlInterval
)
{
    if (debug_)
    {
        Info<< "Write control interval: " << writeControlInterval << endl;
    }

    writeControlInterval_ = writeControlInterval;
}


void AFC::Properties::insertDeltaT
(
    const scalar& deltaT 
)
{
    if (debug_)
    {
        Info<< "Algorithm deltaT: " << deltaT << endl;
    }

    deltaT_ = deltaT;
}


void AFC::Properties::insertPressure
(
    const scalar& pressure
)
{
    p_ = pressure;
}


void AFC::Properties::inputMol()
{
    inputMol_ = true;
}


void AFC::Properties::inputMass()
{
    inputMass_ = true;
}


// * * * * * * * * * * * * * * * Other functions * * * * * * * * * * * * * * //

void AFC::Properties::check()
{
    if (nZPoints_ == 0)
    {
        FatalError
        (
            "    No mixtureFractionPoints defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if (nZvarPoints_ == 0)
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

    if (sDRs_.empty())
    {
        FatalError
        (
            "    No scalar dissipation rates defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if
    (
        oxidizerX_.empty()
     && oxidizerY_.empty()
    )
    {
        FatalError
        (
            "    No oxidizer species defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if
    (
        fuelX_.empty()
     && fuelY_.empty()
    )
    {
        FatalError
        (
            "    No fuel species defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    scalar sum{0};

    scalar epsilon{1e-12};

    if (inputMol_)
    {
        if (debug_)
        {
            Info<< "Input of fraction is mol.\n" << endl;
        }

        //- Oxidizer
        forAll(speciesOxidizer_, s)
        {
            sum += oxidizerX_.at(s);
        } 

        //- Double comparison; with epsilon
        if ((1.0 - sum) > epsilon)
        {
            FatalError
            (
                "    Sum of mol fraction of oxidizer is not 1 ("
                + std::to_string(sum) + ").\n"
                "    Please check your fraction composition for the oxidizer "
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
            FatalError
            (
                "    Sum of mol fraction of fuel is not 1 ("
                + std::to_string(sum) + ").\n"
                "    Please check your fraction composition for the fuel "
                "stream in your afcDict.",
                __FILE__,
                __LINE__
            );
        }
    }
    else if (inputMass_)
    {
        if (debug_)
        {
            Info<< "Input of fraction is mass.\n" << endl;
        }

        //- Oxidizer
        forAll(speciesOxidizer_, s)
        {
            sum += oxidizerY_.at(s);
        } 

        //- Double comparison; with epsilon
        if ((1.0 - sum) > epsilon)
        {
            FatalError
            (
                "    Sum of mass fraction of oxidizer is not 1 ("
                + std::to_string(sum) + ").\n"
                "    Please check your fraction composition for the oxidizer "
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
            FatalError
            (
                "    Sum of mass fraction of fuel is not 1 ("
                + std::to_string(sum) + ").\n"
                "    Please check your fraction composition for the fuel "
                "stream in your afcDict.",
                __FILE__,
                __LINE__
            );
        }
    }
    else
    {
        FatalError
        (
            "    Either mol fraction or mass fraction defined for oxidizer"
            " stream.",
            __FILE__,
            __LINE__
        );
    }

    //- Check if pressure is set
    if (p_ <= 0)
    {
        FatalError
        (
            "    Pressure is not set or set not correct.\n"
            "    Please set the pressure in afcDict.",
            __FILE__,
            __LINE__
        );
    }


    // Algorihm control check

    // Output summary
}


void AFC::Properties::XtoY()
{
    scalar YbyM{0};

    //- get all species from chemistry 
    const wordList& speciesChem = chemistry_.species();

    forAll(speciesChem, speciesI)
    {
        YbyM +=
            thermo_.MW(speciesI)
          * oxidizerX_[speciesI];
    } 

    //- Oxidizer
    {    
        forAll(speciesOxidizer_, speciesI)
        {
            oxidizerY_[speciesI] =
                (oxidizerX_.find(speciesI)->second
              / thermo_.MW(speciesI))
              / YbyM;
        }
    }

    //- Fuel
    {
        forAll(speciesFuel_, speciesI)
        {
            fuelY_[speciesI] =
                (fuelX_.find(speciesI)->second
              / thermo_.MW(speciesI))
              / YbyM;
        }
    }
}


/*void AFC::Properties::YtoX
(
    const word tmp
)
{
    scalar YbyM{0};

    //- get all species from chemistry 
    const wordList& speciesChem = chemistry_.species();

    forAll(speciesChem, s)
    {
        YbyM += thermo_.MW(s) * oxidizerX_.at(s);
    } 

    if (tmp == "O")
    {    
        forAll(speciesOxidizer_, s)
        {
            oxidizerY_[s] = oxidizerX_.at(s) / YbyM;
        }
    }

    if (tmp == "F")
    {
        forAll(speciesFuel_, s)
        {
            fuelY_[s] = fuelX_.at(s) / YbyM;
        }
    }
}*/


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

AFC::wordList AFC::Properties::speciesOxidizer() const
{
    return speciesOxidizer_;
}


AFC::map<AFC::word, AFC::scalar> AFC::Properties::oxidizerCompMol() const
{
    return oxidizerX_;
}


AFC::map<AFC::word, AFC::scalar> AFC::Properties::oxidizerCompMass() const
{
    return oxidizerY_;
}


AFC::wordList AFC::Properties::speciesFuel() const
{
    return speciesFuel_;
}


AFC::map<AFC::word, AFC::scalar> AFC::Properties::fuelCompMol() const
{
    return fuelX_;
}


AFC::map<AFC::word, AFC::scalar> AFC::Properties::fuelCompMass() const
{
    return fuelY_;
}


AFC::scalarField AFC::Properties::defects() const
{
    return defects_;
}


AFC::scalarField AFC::Properties::sDRs() const
{
    return sDRs_;
}


int AFC::Properties::nZPoints() const
{
    return nZPoints_;
}


int AFC::Properties::nZvarPoints() const
{
    return nZvarPoints_;
}


AFC::scalar AFC::Properties::oxidizerTemperature() const
{
    return TOxidizer_;
}


AFC::scalar AFC::Properties::fuelTemperature() const
{
    return TFuel_;
}


AFC::scalar AFC::Properties::runTime() const
{
    return runTime_;
}


AFC::scalar AFC::Properties::deltaT() const
{
    return deltaT_;
}


unsigned int AFC::Properties::nDefects() const
{
    return defects_.size();
}


AFC::scalar AFC::Properties::defect
(
    const int& defectNo
) const
{
    return defects_[defectNo];
}


AFC::scalar AFC::Properties::p() const
{
    return p_;
}


AFC::word AFC::Properties::input() const
{
    if
    (
        !inputMol_
     && !inputMass_
    )
    {
        FatalError
        (
            "    Input problems for mass or mol fraction",
            __FILE__,
            __LINE__
        );
    }

    word ret{"none"};

    if (inputMol_)
    {
        ret = "mol";
    }
    else if (inputMass_)
    {
        ret = "mass";
    }

    return ret;
}


// ************************************************************************* //
