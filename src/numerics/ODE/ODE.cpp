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
#include "ODE.hpp"
#include "stepStatus.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::ODE::ODE
(
    const Chemistry& chem
)
:
    StepStatus(1),

    chem_(chem)
{
    if (debug_)
    {
        Info<< "Constructor ODE \n" << endl;
    }

    solver_ = new Seulex(chem); 
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::ODE::~ODE()
{
    if (debug_)
    {
        Info<< "Destructor ODE \n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

void AFC::ODE::derivative
(
    const scalar dt,
    const map<word, scalar>& c0
)
{
    //- Species      
    const wordList& species = chem_.species();

    forAll(species, s)
    {
        c_[s] = max(scalar(0), c0.at(s));
    }

    forAll(species, s)
    {
        dcdt_[s] = chem_.omega(
    }
}


void AFC::ODE::solve
(
    const scalar T,
    const scalar p,
    map<word, scalar>& c,
    const scalar dt,
    scalar& dtTry
)
{
    size_t iterChem{0};

    scalar timeLeft = dt;

    //- Copy of old concentrations
    const map<word, scalar> c0 = c;

    //- TODO SMALL
    while (timeLeft > 1e-15)
    {
        //- Time step of actual chemistry iteration
        scalar dt = 0.00001;

        //- Calculate the derivative dcdt
        derivative(dt, c0);

        timeLeft -= dt;

        iterChem++;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
