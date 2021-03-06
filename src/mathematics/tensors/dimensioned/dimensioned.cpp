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

#include "dimensioned.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<typename Type>
TKC::Dimensioned<Type>::Dimensioned
(
    const size_t kg,
    const size_t m,
    const size_t s,
    const size_t K
)
:
    value_(4)
{}


template<typename Type>
TKC::Dimensioned<Type>::Dimensioned
(
    const scalarField& dimension 
)
: 
    value_(4)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<typename Type>
TKC::Dimensioned<Type>::~Dimensioned()
{}


// * * * * * * * * * * * * * * Operator Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * //


// ************************************************************************* //
