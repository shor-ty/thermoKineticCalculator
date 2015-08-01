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
{
    Info<< "Destructor TransportData\n";
}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //



// * * * * * * * * * Insert functions from TransportReader:: * * * * * * * * //

void AFC::TransportData::insertSpecies
(
    const word& species
)
{
    species_.push_back(species);
}


void AFC::TransportData::insertGeoConfig
(
    const int& geoConfig
)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == geoConfig_.size())
    {
        geoConfig_[species_[species_.size()-1]] = geoConfig;
    }
    else
    {
        FatalError
        (
            "    geoConfig_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::TransportData::insertLenJonPot
(
    const scalar& lenJonPot
)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == lenJonPot_.size())
    {
        lenJonPot_[species_[species_.size()-1]] = lenJonPot;
    }
    else
    {
        FatalError
        (
            "    lenJonPot_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::TransportData::insertLenJonCollDia
(
    const scalar& lenJonCollDia
)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == lenJonCollDia_.size())
    {
        lenJonCollDia_[species_[species_.size()-1]] = lenJonCollDia;
    }
    else
    {
        FatalError
        (
            "    lenJonCollDia_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::TransportData::insertDipMom
(
    const scalar& dipMom
)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == dipMom_.size())
    {
        dipMom_[species_[species_.size()-1]] = dipMom;
    }
    else
    {
        FatalError
        (
            "    dipMom_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::TransportData::insertPol
(
    const scalar& pol
)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == pol_.size())
    {
        pol_[species_[species_.size()-1]] = pol;
    }
    else
    {
        FatalError
        (
            "    pol_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::TransportData::insertRotRelCollNumb
(
    const scalar& rotRelCollNumb
)
{
    //- Species_ list must have one element more in this list 
    if (species_.size()-1 == rotRelCollNumb_.size())
    {
        rotRelCollNumb_[species_[species_.size()-1]] = rotRelCollNumb;
    }
    else
    {
        FatalError
        (
            "    rotRelCollNumb_.size() is not equal to species_.size()-1.\n"
            "    Some error occur.",
            __FILE__,
            __LINE__
        );
    }
}


// * * * * * * * * * * * * * Setter bool functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Return functions  * * * * * * * * * * * * * * //

AFC::wordList AFC::TransportData::species() const
{
    return species_;
}


// ************************************************************************* //
