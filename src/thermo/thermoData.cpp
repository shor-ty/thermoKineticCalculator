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

#include "thermoData.hpp"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

AFC::ThermoData::ThermoData
(
    const bool& thermo
)
:
    thermo_{thermo}

{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

AFC::ThermoData::~ThermoData()
{
    Info<< "Destructor ThermoData\n";
}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //



// * * * * * * * * * Insert functions from ThermoReader:: * * * * * * * * //

void AFC::ThermoData::insertSpecies
(
    const word& species
)
{
    species_.push_back(species);
}


void AFC::ThermoData::insertChemicalFormula
(
    const word& chemFormula 
)
{
    formula_.push_back(chemFormula);
}


void AFC::ThermoData::insertAtomAndFactor
(
    const word& atom,
    const unsigned int& factor
)
{
    //- Using normal names
    const word& actualSpecies = species_[species_.size()-1];

    //- Insert atom to correct position
    speciesAtoms_[actualSpecies].push_back(atom);

    //- Insert multiplication factor
    atomFactors_[actualSpecies].push_back(factor);
}


void AFC::ThermoData::insertMolecularWeight
(
    const scalar& MW 
)
{
    //- Species_ list must have one element more in list than MW_
    if (species_.size()-1 == MW_.size())
    {
        MW_[species_[species_.size()-1]] = MW;
    }
    else
    {
        FatalError
        (
            "    MW_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::ThermoData::insertPhase
(
    const word& phase
)
{
    
    //- Species_ list must have one element more in list than phase_
    if (species_.size()-1 == phase_.size())
    {
        phase_[species_[species_.size()-1]] = phase;
    }
    else
    {
        FatalError
        (
            "    phase_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::ThermoData::insertLT
(
    const scalar& LT
)
{
    
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == LT_.size())
    {
        LT_[species_[species_.size()-1]] = LT;
    }
    else
    {
        FatalError
        (
            "    LT_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::ThermoData::insertHT
(
    const scalar& HT
)
{
    
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == HT_.size())
    {
        HT_[species_[species_.size()-1]] = HT;
    }
    else
    {
        FatalError
        (
            "    HT_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::ThermoData::insertCT
(
    const scalar& CT
)
{
    
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == CT_.size())
    {
        CT_[species_[species_.size()-1]] = CT;
    }
    else
    {
        FatalError
        (
            "    CT_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::ThermoData::insertNASACoeffsHT
(
    const scalar& pc
)
{
    //- Insert value
    NASACoeffsHT_[species_[species_.size()-1]].push_back(pc);
}


void AFC::ThermoData::insertNASACoeffsLT
(
    const scalar& pc
)
{
    //- Insert value
    NASACoeffsLT_[species_[species_.size()-1]].push_back(pc);
}


// * * * * * * * * * * * Insert functions from Thermo::  * * * * * * * * * * //

void AFC::ThermoData::p
(
    const scalar& pressure
) 
{
    p_ = pressure;
}


// * * * * * * * * * * * * * Setter bool functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Return functions  * * * * * * * * * * * * * * //

AFC::wordList AFC::ThermoData::species() const
{
    return species_;
}


AFC::scalar AFC::ThermoData::MW
(
    const word& species
) const
{
//    return MW_.find(species)->second;
    return MW_.at(species);
}


AFC::scalar AFC::ThermoData::LT
(
    const word& species
) const
{
    return LT_.at(species);
}


AFC::scalar AFC::ThermoData::CT
(
    const word& species
) const
{
    return CT_.at(species);
}


AFC::scalar AFC::ThermoData::HT
(
    const word& species
) const
{
    return HT_.at(species);
}


AFC::scalarField AFC::ThermoData::NASACoeffsHT
(
    const word& species
) const
{
    return NASACoeffsHT_.at(species); 
}


AFC::scalarField AFC::ThermoData::NASACoeffsLT
(
    const word& species
) const
{
    return NASACoeffsLT_.at(species); 
}


AFC::scalar AFC::ThermoData::p() const
{
    return p_;
}


AFC::wordList AFC::ThermoData::speciesAtoms
(
    const word& species
) const
{
    return speciesAtoms_.at(species);
}


AFC::scalarList AFC::ThermoData::atomFactors
(
    const word& species
) const
{
    return atomFactors_.at(species);
}


// ************************************************************************* //
