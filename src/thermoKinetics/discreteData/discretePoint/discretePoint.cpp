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

#include "discretePoint.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::DiscretePoint::DiscretePoint()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::DiscretePoint::~DiscretePoint()
{}


// * * * * * * * * * * * * * * * Insert Functions  * * * * * * * * * * * * * //

void TKC::DiscretePoint::C(const word species, const scalar value)
{
    C_[species] = value;
}


void TKC::DiscretePoint::Y(const word species, const scalar value)
{
    Y_[species] = value;
}


void TKC::DiscretePoint::X(const word species, const scalar value)
{
    X_[species] = value;
}


void TKC::DiscretePoint::T(const scalar value)
{
    T_ = value;
}


void TKC::DiscretePoint::rho(const scalar value)
{
    rho_ = value;
}


void TKC::DiscretePoint::MMW(const scalar value)
{
    MMW_ = value;
}


// * * * * * * * * * * * * * * * Return Funcitons  * * * * * * * * * * * * * //

const TKC::map<TKC::word, TKC::scalar>& TKC::DiscretePoint::C() const
{
    return C_;
}


const TKC::scalar& TKC::DiscretePoint::C(const word species) const
{
    return C_.at(species);
}


const TKC::map<TKC::word, TKC::scalar>& TKC::DiscretePoint::Y() const
{
    return Y_;
}


const TKC::scalar& TKC::DiscretePoint::Y(const word species) const
{
    return Y_.at(species);
}


const TKC::map<TKC::word, TKC::scalar>& TKC::DiscretePoint::X() const
{
    return X_;
}


const TKC::scalar& TKC::DiscretePoint::X(const word species) const
{
    return X_.at(species);
}


const TKC::scalar& TKC::DiscretePoint::T() const
{
    return T_;
}


const TKC::scalar& TKC::DiscretePoint::rho() const
{
    return rho_;
}


const TKC::scalar& TKC::DiscretePoint::MMW() const
{
    return MMW_;
}


// ************************************************************************* //
