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

#include "chemistry.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Chemistry::Chemistry()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Chemistry::~Chemistry()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AFC::Chemistry::chemistryReader
(
    const string& fileName 
)
{
    pCR_ = smartPtr<ChemistryReader>(new ChemistryReader(fileName, this));
    pCD_ = smartPtr<ChemistryData>(new ChemistryData);
}


void AFC::Chemistry::readChemistry()
{
    pCR_->readChemistry();
}


// * * * * * * * * * Insert functions for ChemistryReader::  * * * * * * * * //

void AFC::Chemistry::insertElements
(
    const word& element
)
{
    pCD_->insertElements(element);
}


void AFC::Chemistry::insertSpecies
(
    const word& species
)
{
    pCD_->insertSpecies(species);
}


void AFC::Chemistry::insertElementarReaction
(
    const string& reaction
)
{
    pCD_->insertElementarReaction(reaction);
}


void AFC::Chemistry::insertArrheniusCoeffs
(
    const scalar& coeff_1,
    const scalar& coeff_2,
    const scalar& coeff_3
)
{
    pCD_->insertArrheniusCoeffs
        (
            coeff_1,
            coeff_2,
            coeff_3
        );
}


void AFC::Chemistry::insertLOWCoeffs
(
    const scalar& coeff,
    const unsigned int& coeffNo
)
{
    pCD_->insertLOWCoeffs(coeff, coeffNo);
}


void AFC::Chemistry::insertTROECoeffs
(
    const scalar& coeff,
    const unsigned int& coeffNo
)
{
    pCD_->insertTROECoeffs(coeff, coeffNo);
}


void AFC::Chemistry::insertSRICoeffs
(
    const scalar& coeff,
    const unsigned int& coeffNo
)
{
    pCD_->insertSRICoeffs(coeff, coeffNo);
}


void AFC::Chemistry::insertMvalue
(
    const scalar& value
)
{
    pCD_->insertMvalue(value);
}


void AFC::Chemistry::insertMcomp
(
    const string& comp
)
{
    pCD_->insertMcomp(comp);
}


void AFC::Chemistry::incrementDuplicated()
{
    pCD_->incrementDuplicated();
}


void AFC::Chemistry::incrementReac()
{
    pCD_->incrementReac();
}


void AFC::Chemistry::incrementMatrixesVectors()
{
    pCD_->incrementMatrixesVectors();
}


// * * * * * * * * * * * * * Setter bool functions * * * * * * * * * * * * * //

void AFC::Chemistry::setKB()
{
    pCD_->setKB();
}


void AFC::Chemistry::setTBR()
{
    pCD_->setTBR();
}


void AFC::Chemistry::setLOW()
{
    pCD_->setLOW();
}


void AFC::Chemistry::setTROE()
{
    pCD_->setTROE();
}


void AFC::Chemistry::setSRI()
{
    pCD_->setSRI();
}


void AFC::Chemistry::setENHANCE()
{
    pCD_->setENHANCE();
}


// ************************************************************************* //
