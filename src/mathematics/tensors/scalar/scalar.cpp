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

#include "definitions.hpp"
#include "scalar.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::Scalar::Scalar()
:
    Tensor(0, 0)
{}


/*TKC::Scalar::Scalar
(
    const Scalar& vec
)
:
    Tensor(vec.rows(), vec.cols(), vec.values())
{}


TKC::Scalar::Scalar
(
    const size_t row,
    const scalar value
)
:
    Tensor(row, 1, value)
{}


TKC::Scalar::Scalar
(
    const scalarField& sF
)
:
    Tensor(sF.size(), 1, sF)
{}*/


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::Scalar::~Scalar()
{}


// * * * * * * * * * * * * * * Operator Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * //

//TKC::Scalar TKC::Scalar

// ************************************************************************* //
