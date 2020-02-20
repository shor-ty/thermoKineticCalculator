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

#include "transportData.hpp"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

AFC::TransportData::TransportData()
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

AFC::TransportData::~TransportData()
{}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * Fitting polynomials coefficients  * * * * * * * * * //

void AFC::TransportData::viscosityPolyCoeffs
(
    const word species,
    const Vector& x
)
{
    //- Solution vector x containts polynomial coefficients D, C, B and A
    viscosity_[species] = x.values();
}


AFC::scalarField
AFC::TransportData::viscosityPolyCoeffs(const word species) const
{
    return viscosity_.at(species);
}


void AFC::TransportData::thermalConductivityPolyCoeffs
(
    const word species,
    const Vector& x
)
{
    //- Solution vector x containts polynomial coefficients D, C, B and A
    thermalConductivity_[species] = x.values();
}


AFC::scalarField
AFC::TransportData::thermalConductivityPolyCoeffs(const word species) const
{
    return thermalConductivity_.at(species);
}


void AFC::TransportData::binaryDiffusivityPolyCoeffs
(
    const word species1,
    const word species2,
    const Vector& x
)
{
    //- Check which species we have (ID)
    const wordList& chemistrySpecies_ = chemistrySpecies();

    size_t ID{0};

    forAll(chemistrySpecies_, s)
    {
        if (s == species1)
        {
            break;
        }

        ++ID;
    }

    //- If size equal to ID -> is available and we just add the
    //  species2 and polynomial coefficents to the map
    if (binaryDiffusivity_.size() - 1 == ID)
    {
        binaryDiffusivity_[ID][species2] = x.values();
    }
    //- If ID == size()+1; not available and hence we first have to
    //  extend the List by one
    else if (binaryDiffusivity_.size() == ID)
    {
        binaryDiffusivity_.push_back(map<word, scalarField>());

        binaryDiffusivity_[ID][species2] = x.values();
    }
    else
    {
        ErrorMsg
        (
            "    Building the binaryDiffusivity_ field went wrong\n"
            "    Please check the code or ask for help",
            __FILE__,
            __LINE__
        );
    }
}


AFC::scalarField AFC::TransportData::binaryDiffusivityPolyCoeffs
(
    const word species1,
    const word species2
) const
{
    //- First find the ID of first species
    const wordList& chemistrySpecies_ = chemistrySpecies();

    size_t ID{0};

    forAll(chemistrySpecies_, s)
    {
        if (s == species1)
        {
            break;
        }

        ++ID;
    }
    
    //- Return the coefficients
    return binaryDiffusivity_[ID].at(species2);
}


// * * * * * * * * * Insert functions from TransportReader:: * * * * * * * * //

void AFC::TransportData::insertSpecies(const word species)
{
    species_.push_back(species);
}


void AFC::TransportData::insertChemistrySpecies
(
    const wordList& chemistrySpecies
)
{
    chemistrySpecies_ = chemistrySpecies;
}


void AFC::TransportData::insertGeoConfig(const int geoConfig)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == geoConfig_.size())
    {
        geoConfig_[species_[species_.size()-1]] = geoConfig;
    }
    else
    {
        ErrorMsg
        (
            "    geoConfig_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::TransportData::insertLenJonPot(const scalar lenJonPot)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == lenJonPot_.size())
    {
        lenJonPot_[species_[species_.size()-1]] = lenJonPot;
    }
    else
    {
        ErrorMsg
        (
            "    lenJonPot_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::TransportData::insertLenJonCollDia(const scalar lenJonCollDia)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == lenJonCollDia_.size())
    {
        lenJonCollDia_[species_[species_.size()-1]] = lenJonCollDia;
    }
    else
    {
        ErrorMsg
        (
            "    lenJonCollDia_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::TransportData::insertDipMom(const scalar dipMom)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == dipMom_.size())
    {
        dipMom_[species_[species_.size()-1]] = dipMom;
    }
    else
    {
        ErrorMsg
        (
            "    dipMom_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::TransportData::insertAlpha(const scalar alpha)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == alpha_.size())
    {
        alpha_[species_[species_.size()-1]] = alpha;
    }
    else
    {
        ErrorMsg
        (
            "    alpha_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::TransportData::insertZRot298(const scalar ZRot298)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == ZRot298_.size())
    {
        ZRot298_[species_[species_.size()-1]] = ZRot298;
    }
    else
    {
        ErrorMsg
        (
            "    ZRot298_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::TransportData::insertBinarySpeciesCombinations
(
    const word parentSpecies,
    const word childSpecies 
)
{
    binarySpeciesCombinations_[parentSpecies].push_back(childSpecies);
}


// * * * * * * * * * * * * Insert functions from afc.cpp * * * * * * * * * * //

void AFC::TransportData::chemicalFormula
(
    const wordList& chemicalFormula
)
{
    chemicalFormula_ = chemicalFormula;
}


// * * * * * * * * * * * Update and Manipulation Functions * * * * * * * * * //

void AFC::TransportData::binarySpeciesCombinations()
{
    const wordList& species_ = species();

    bool found{false};

    forAll(species_, firstSpecies)
    {
        forAll(species_, secondSpecies)
        {
            if (found)
            {
                insertBinarySpeciesCombinations(firstSpecies, secondSpecies);
            }

            if (firstSpecies == secondSpecies)
            {
                found = true;
            }
        }

        found = false;
    }
}


// * * * * * * * * * * * * * * Return functions  * * * * * * * * * * * * * * //

AFC::wordList AFC::TransportData::species() const
{
    return species_;
}


AFC::wordList AFC::TransportData::chemicalFormula() const
{
    return chemicalFormula_;
}


AFC::word AFC::TransportData::chemicalFormula(const word species) const
{
    //- TODO use map to speed up 
    int ID{0};
    
    //- Search id
    forAll(species_, s)
    {
        if (s == species)
        {
            break;
        }

        ++ID;
    }

    return chemicalFormula_[ID];
}


AFC::wordList AFC::TransportData::chemistrySpecies() const
{
    return chemistrySpecies_;
}


int AFC::TransportData::geometricalConfig(const word species) const
{
    return geoConfig_.at(species);
}


AFC::scalar AFC::TransportData::LJCD(const word species) const
{
    return lenJonCollDia_.at(species);
}


AFC::scalar AFC::TransportData::LJP(const word species) const
{
    return lenJonPot_.at(species);
}


AFC::scalar AFC::TransportData::muk(const word species) const
{
    return dipMom_.at(species);
}


AFC::scalar AFC::TransportData::alpha(const word species) const
{
    return alpha_.at(species);
}


AFC::scalar AFC::TransportData::ZRot298(const word species) const
{
    return ZRot298_.at(species);
}


// ************************************************************************* //
