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

#include "thermoData.hpp"
#include "thermoReader.hpp"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

TKC::ThermoData::ThermoData
(
    const string fileName,
    const bool thermoInChemistry
)
:
    thermoInChemistry_(thermoInChemistry),
    // TODO
    p_(12)
{
    if (debug_)
    {
        Info<< "ThermoData Constructor\n" << endl;
    }

    //- Read the thermodynamic data
    {
        ThermoReader reader(fileName, *this);

        reader.read();
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

TKC::ThermoData::~ThermoData()
{
    if (debug_)
    {
        Info<< "ThermoData Destructor\n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //



// * * * * * * * * * Insert functions from ThermoReader:: * * * * * * * * //

void TKC::ThermoData::setSpecies(const word species)
{
    species_.push_back(species);
}


void TKC::ThermoData::setChemicalFormula(const word chemFormula)
{
    formula_.push_back(chemFormula);
}


void TKC::ThermoData::setElementAndAtoms
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
    elementAtoms_[actualSpecies].push_back(factor);
}


void TKC::ThermoData::setMolecularWeight(const scalar MW) 
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


void TKC::ThermoData::setPhase(const word phase)
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


void TKC::ThermoData::setLT(const scalar LT)
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


void TKC::ThermoData::setHT(const scalar HT)
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


void TKC::ThermoData::setCT(const scalar CT)
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


void TKC::ThermoData::setNASACoeffsHT(const scalar coeff)
{
    //- Insert value
    NASACoeffsHT_[species_[species_.size()-1]].push_back(coeff);
}


void TKC::ThermoData::setNASACoeffsLT(const scalar coeff)
{
    //- Insert value
    NASACoeffsLT_[species_[species_.size()-1]].push_back(coeff);
}


void TKC::ThermoData::updateElementsAndFactors()
{
    //- Temporary map that contains all elements and factors of
    //  the actual species
    map<word, scalar> tmp;

    //- Actual species
    const word species = species_[species_.size()-1];

    const wordList& elements = elementsInSpecies(species);
    const scalarList& factors = elementAtoms(species);

    forEach(elements, a)
    {
        tmp[elements[a]] = factors[a];
    }

    elements_[species] = tmp;
}


// * * * * * * * * * * * Insert functions from Thermo::  * * * * * * * * * * //

void TKC::ThermoData::p(const scalar pressure) 
{
    p_ = pressure;
}


// * * * * * * * * * * * * * Setter bool functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Return functions  * * * * * * * * * * * * * * //

TKC::scalar TKC::ThermoData::p() const
{
    return p_;
}


const TKC::wordList TKC::ThermoData::species() const
{
    return species_;
}


const TKC::wordList TKC::ThermoData::formula() const
{
    return formula_;
}


const TKC::wordList
TKC::ThermoData::elementsInSpecies(const word species) const
{
    return elementsInSpecies_.at(species);
}


const TKC::wordList
TKC::ThermoData::elementsInSpeciesChem(const word species) const
{
    NotImplemented(__FILE__, __LINE__);

    return wordList(0);
}


const TKC::scalarList TKC::ThermoData::elementAtoms(const word species) const
{
    return elementAtoms_.at(species);
}


const TKC::map<TKC::word, TKC::scalar>
TKC::ThermoData::elementAtomsMap(const word species) const
{
    return elements_.at(species);
}


const TKC::map<TKC::word, TKC::scalar>
TKC::ThermoData::elementAtomsChem(const word species) const
{
    NotImplemented(__FILE__, __LINE__);

    map<word, scalar> temp;
    temp["NONE"] = scalar(0);

    return temp;
}


const TKC::map<TKC::word, TKC::scalar> TKC::ThermoData::MW() const
{
    return MW_;
}


TKC::scalar TKC::ThermoData::MW(const word species) const
{
    return MW_.at(species);
}


const TKC::map<TKC::word, TKC::word> TKC::ThermoData::phase() const
{
    return phase_;
}


const TKC::word TKC::ThermoData::phase(const word species) const
{
    return phase_.at(species);
}


TKC::scalar TKC::ThermoData::LT(const word species) const
{
    return LT_.at(species);
}


TKC::scalar TKC::ThermoData::CT(const word species) const
{
    return CT_.at(species);
}


TKC::scalar TKC::ThermoData::HT(const word species) const
{
    return HT_.at(species);
}


const TKC::scalarField TKC::ThermoData::NASACoeffsLT(const word species) const
{
    return NASACoeffsLT_.at(species); 
}


const TKC::scalarField TKC::ThermoData::NASACoeffsHT(const word species) const
{
    return NASACoeffsHT_.at(species); 
}



// ************************************************************************* //
