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

#include "typedef.hpp"
#include "scalar.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Scalar::Scalar()
:
    Tensor(0, 0)
{
    if (debug_)
    {
        Info<< "Constructor Scalar() \n" << endl;
    }
}


/*AFC::Scalar::Scalar
(
    const Scalar& vec
)
:
    Tensor(vec.rows(), vec.cols(), vec.values())
{
    if (debug_)
    {
        Info<< "Constructor Scalar(const Scalar&) Copy Const\n" << endl;
    }
}


AFC::Scalar::Scalar
(
    const size_t row,
    const scalar value
)
:
    Tensor(row, 1, value)
{
    if (debug_)
    {
        Info<< "Constructor Scalar (row, value)\n" << endl;
    }
}


AFC::Scalar::Scalar
(
    const scalarField& sF
)
:
    Tensor(sF.size(), 1, sF)
{
    if (debug_)
    {
        Info<< "Constructor Scalar (row, scalarField)\n" << endl;
    }
}*/


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Scalar::~Scalar()
{
    if (debug_)
    {
        Info<< "Destructor Scalar\n" << endl;
    }
}


// * * * * * * * * * * * * * * Operator Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * //

//AFC::Scalar AFC::Scalar

// ************************************************************************* //
