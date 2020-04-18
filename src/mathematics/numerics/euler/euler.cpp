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

#include "definitions.hpp"
#include "euler.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Euler::Euler(const size_t nSpecies)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Euler::~Euler()
{}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

AFC::scalar AFC::Euler::solve
(
    const scalar t0,
    const scalar dt,
    const map<word, scalar>& c0,
    const map<word, scalar>& dcdt,
    map<word, scalar>& c
) const
{
    //- TODO make everything consitent
    //  Map<word, scalar> replace by vector or find some nice solution
    size_t i{0};

    //- Calculate error estimated by the change in the state
    //  and update the state
    forAll(c0, s)
    {
        err_[s.first] = dt * dcdt.at(s.first);
        c[s.first] = err_.at(s.first) * s.second;
    }

    //- Normalize error
    //return normalizeError(c0, c, err_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
