/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
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

AFC::ThermoData::ThermoData(const bool thermo)
:
    thermo_{thermo}

{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

AFC::ThermoData::~ThermoData()
{}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //



// * * * * * * * * * Insert functions from ThermoReader:: * * * * * * * * //

void AFC::ThermoData::insertSpecies(const word species)
{
    species_.push_back(species);
}


void AFC::ThermoData::insertChemicalFormula(const word chemFormula)
{
    formula_.push_back(chemFormula);
}


void AFC::ThermoData::insertElementAndFactor
(
    const word atom,
    const unsigned int factor
)
{
    //- Using normal names
    const word& actualSpecies = species_[species_.size()-1];

    //- Insert atom to correct position
    elementsInSpecies_[actualSpecies].push_back(atom);

    //- Insert multiplication factor
    elementFactors_[actualSpecies].push_back(factor);
}


void AFC::ThermoData::insertMolecularWeight(const scalar MW) 
{
    //- Species_ list must have one element more in list than MW_
    if (species_.size()-1 == MW_.size())
    {
        MW_[species_[species_.size()-1]] = MW;
    }
    else
    {
        ErrorMsg
        (
            "    MW_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::ThermoData::insertPhase(const word phase)
{
    //- Species_ list must have one element more in list than phase_
    if (species_.size()-1 == phase_.size())
    {
        phase_[species_[species_.size()-1]] = phase;
    }
    else
    {
        ErrorMsg
        (
            "    phase_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::ThermoData::insertLT(const scalar LT)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == LT_.size())
    {
        LT_[species_[species_.size()-1]] = LT;
    }
    else
    {
        ErrorMsg
        (
            "    LT_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::ThermoData::insertHT(const scalar HT)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == HT_.size())
    {
        HT_[species_[species_.size()-1]] = HT;
    }
    else
    {
        ErrorMsg
        (
            "    HT_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::ThermoData::insertCT(const scalar CT)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == CT_.size())
    {
        CT_[species_[species_.size()-1]] = CT;
    }
    else
    {
        ErrorMsg
        (
            "    CT_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::ThermoData::insertNASACoeffsHT(const scalar coeff)
{
    //- Insert value
    NASACoeffsHT_[species_[species_.size()-1]].push_back(coeff);
}


void AFC::ThermoData::insertNASACoeffsLT(const scalar coeff)
{
    //- Insert value
    NASACoeffsLT_[species_[species_.size()-1]].push_back(coeff);
}


void AFC::ThermoData::updateElementsAndFactors()
{
    //- Temporary map that contains all elements and factors of
    //  the actual species
    map<word, scalar> tmp;

    //- Actual species
    const word species = species_[species_.size()-1];

    const wordList& elements = elementsInSpecies(species);
    const scalarList& factors = elementFactors(species);

    forEach(elements, a)
    {
        tmp[elements[a]] = factors[a];
    }

    elements_[species] = tmp;
}


// * * * * * * * * * * * Insert functions from Thermo::  * * * * * * * * * * //

void AFC::ThermoData::p(const scalar pressure) 
{
    p_ = pressure;
}


// * * * * * * * * * * * * * Setter bool functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Return functions  * * * * * * * * * * * * * * //

AFC::scalar AFC::ThermoData::p() const
{
    return p_;
}


const AFC::wordList AFC::ThermoData::species() const
{
    return species_;
}


const AFC::wordList AFC::ThermoData::formula() const
{
    return formula_;
}


const AFC::wordList
AFC::ThermoData::elementsInSpecies(const word species) const
{
    return elementsInSpecies_.at(species);
}


const AFC::wordList
AFC::ThermoData::elementsInSpeciesChem(const word species) const
{
    NotImplemented(__FILE__, __LINE__);

    return wordList(0);
}


const AFC::scalarList AFC::ThermoData::elementFactors(const word species) const
{
    return elementFactors_.at(species);
}


const AFC::map<AFC::word, AFC::scalar>
AFC::ThermoData::elementFactorsMap(const word species) const
{
    return elements_.at(species);
}


const AFC::map<AFC::word, AFC::scalar>
AFC::ThermoData::elementFactorsChem(const word species) const
{
    NotImplemented(__FILE__, __LINE__);

    map<word, scalar> temp;
    temp["NONE"] = scalar(0);

    return temp;
}


const AFC::map<AFC::word, AFC::scalar> AFC::ThermoData::MW() const
{
    return MW_;
}


AFC::scalar AFC::ThermoData::MW(const word species) const
{
    return MW_.at(species);
}


const AFC::map<AFC::word, AFC::word> AFC::ThermoData::phase() const
{
    return phase_;
}


const AFC::word AFC::ThermoData::phase(const word species) const
{
    return phase_.at(species);
}


AFC::scalar AFC::ThermoData::LT(const word species) const
{
    return LT_.at(species);
}


AFC::scalar AFC::ThermoData::CT(const word species) const
{
    return CT_.at(species);
}


AFC::scalar AFC::ThermoData::HT(const word species) const
{
    return HT_.at(species);
}


const AFC::scalarField AFC::ThermoData::NASACoeffsLT(const word species) const
{
    return NASACoeffsLT_.at(species); 
}


const AFC::scalarField AFC::ThermoData::NASACoeffsHT(const word species) const
{
    return NASACoeffsHT_.at(species); 
}


// ************************************************************************* //
