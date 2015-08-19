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
    const Thermo& thermo,
    const Chemistry& chemistry
)
:
    thermo_(thermo),

    chemistry_(chemistry)
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
    const scalar& molFraction,
    const bool& lastEntry
)
{
    if (debug)
    {
        Info<< "Oxidizer species: " << species << "  " << molFraction << endl;
    }

    speciesOxidizer_.push_back(species);

    oxidizerX_[species] = molFraction;

    if (lastEntry)
    {
        XtoY("O");
    }
}


void AFC::Properties::insertCompositionFuelMol
(
    const word& species,
    const scalar& molFraction,
    const bool& lastEntry
)
{
    if (debug)
    {
        Info<< "Fuel species: " << species << "  " << molFraction << endl;
    }

    speciesFuel_.push_back(species);

    fuelX_[species] = molFraction;

    if (lastEntry)
    {
        XtoY("F");
    }
}


void AFC::Properties::insertRunTime
(
    const scalar& runTime 
)
{
    if (debug)
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
    if (debug)
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
    if (debug)
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
    if (debug)
    {
        Info<< "Algorithm deltaT: " << deltaT << endl;
    }

    deltaT_ = deltaT;
}


// * * * * * * * * * * * * * * * Other functions * * * * * * * * * * * * * * //

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

    if (sDRs_.empty())
    {
        FatalError
        (
            "    No scalar dissipation rates defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if (oxidizerX_.empty())
    {
        FatalError
        (
            "    No oxidizer species defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if (fuelX_.empty())
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
        sum += oxidizerX_[speciesOxidizer_[i]];
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
        sum += fuelX_[speciesFuel_[i]];
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

    // Algorihm control check

        //- TODO
}


void AFC::Properties::XtoY
(
    const word tmp
)
{
    scalar YbyM{0};

    //- get all species from chemistry 
    const wordList& speciesChem = chemistry_.species();

    forAll(speciesChem, speciesI)
    {
        YbyM +=
            thermo_.MW(speciesChem[speciesI])
          * oxidizerX_[speciesChem[speciesI]];
    } 

    if (tmp == "O")
    {    
        forAll(speciesOxidizer_, speciesI)
        {
            oxidizerY_[speciesOxidizer_[speciesI]] =
                (oxidizerX_.find(speciesOxidizer_[speciesI])->second
              / thermo_.MW(speciesOxidizer_[speciesI]))
              / YbyM;
        }
    }

    if (tmp == "F")
    {
        forAll(speciesFuel_, speciesI)
        {
            fuelY_[speciesFuel_[speciesI]] =
                (fuelX_.find(speciesFuel_[speciesI])->second
              / thermo_.MW(speciesFuel_[speciesI]))
              / YbyM;
        }
    }
}


void AFC::Properties::YtoX
(
    const word tmp
)
{
    scalar YbyM{0};

    //- get all species from chemistry 
    const wordList& speciesChem = chemistry_.species();

    forAll(speciesChem, speciesI)
    {
        const word& speciesI_ = speciesChem[speciesI];

        YbyM += thermo_.MW(speciesI_) * oxidizerX_[speciesChem[speciesI]];
    } 

    if (tmp == "O")
    {    
        forAll(speciesOxidizer_, speciesI)
        {
            oxidizerY_[speciesOxidizer_[speciesI]] =
                oxidizerX_.find(speciesOxidizer_[speciesI])->second
              / YbyM;
        }
    }

    if (tmp == "F")
    {
        forAll(speciesFuel_, speciesI)
        {
            fuelY_[speciesFuel_[speciesI]] =
                fuelX_.find(speciesFuel_[speciesI])->second
              / YbyM;
        }
    }
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

AFC::wordList AFC::Properties::speciesOxidizer() const
{
    return speciesOxidizer_;
}


AFC::map<AFC::word, AFC::scalar> AFC::Properties::oxidizerCompMol() const
{
    return oxidizerX_;
}


AFC::wordList AFC::Properties::speciesFuel() const
{
    return speciesFuel_;
}


AFC::map<AFC::word, AFC::scalar> AFC::Properties::fuelCompMol() const
{
    return fuelX_;
}


AFC::scalarField AFC::Properties::defects() const
{
    return defects_;
}


AFC::scalarField AFC::Properties::sDRs() const
{
    return sDRs_;
}


int AFC::Properties::mfPoints() const
{
    return mfPoints_;
}


int AFC::Properties::vmfPoints() const
{
    return vmfPoints_;
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


// ************************************************************************* //
