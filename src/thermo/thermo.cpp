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

#include "thermo.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Thermo::Thermo
(
    const string& fileName,
    const bool& thermo
)
:
    thermoData_(thermo)
{
    ThermoReader thermoReader(fileName);

    thermoReader.read(thermoData_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Thermo::~Thermo()
{}


// * * * * * * * * * * * * * * * Insert Functions  * * * * * * * * * * * * * //

void AFC::Thermo::p
(
    const scalar& pressure
)
{
    thermoData_.p(pressure);
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

AFC::wordList AFC::Thermo::species() const
{
    return thermoData_.species();
}


AFC::scalar AFC::Thermo::MW
(
    const word& species
) const
{
    return thermoData_.MW(species);
}


AFC::map<AFC::word, AFC::scalar> AFC::Thermo::MW() const
{
    return thermoData_.MW();
}


AFC::scalar AFC::Thermo::MmeanX
(
    const map<word, scalar>& X
) const
{
    return thermoCalc_.MmeanX(X, MW());    
}


AFC::scalar AFC::Thermo::cp
(
    const word& species,
    const scalar& T
) const
{
    return thermoCalc_.cp(species, T, thermoData_);
}


AFC::scalar AFC::Thermo::H
(
    const word& species,
    const scalar& T
) const
{
    return thermoCalc_.H(species, T, thermoData_);
}


AFC::scalar AFC::Thermo::S
(
    const word& species,
    const scalar& T
) const
{
    return thermoCalc_.S(species, T, thermoData_);
}


AFC::scalar AFC::Thermo::G
(
    const word& species,
    const scalar& T
) const
{
    return thermoCalc_.G(species, T, thermoData_);
}


AFC::scalar AFC::Thermo::G
(
    const scalar& H,
    const scalar& S,
    const scalar& T
) const
{
    return thermoCalc_.G(H, S, T);
}


AFC::scalar AFC::Thermo::p() const
{
    return thermoData_.p();
}


AFC::scalar AFC::Thermo::Hf
(
    const word& species,
    const scalar& T
) const
{
    return thermoCalc_.Hf(species, T, thermoData_);
}


AFC::scalar AFC::Thermo::C
(
    const scalar& T
) const
{
    return thermoCalc_.C(p(), T);
}


// * * * * * * * * * * * * * * Summary function  * * * * * * * * * * * * * * //

void AFC::Thermo::summary
(
    ostream& data
) const
{


}


// ************************************************************************* //
