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

#include "typedef.hpp"
#include "rosenbrock.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Rosenbrock::Rosenbrock
(
    Chemistry& chemistry_
)
:
    ODE(chemistry_)
{
    if (debug_)
    {
        Info<< "Constructor for Rosenbrock\n" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Rosenbrock::~Rosenbrock()
{
    if (debug_)
    {
        Info<< "Destructor for Rosenbrock\n" << endl;
    }
}

// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

void AFC::Rosenbrock::solve
(
    const scalar Ta,
    const wordList& species,
    const map<word, scalar>& c1
)
{

    scalar T = 1000;
    map<word, scalar> c = c1;

    Info<< std::setw(10) << "Iter.";
    forAll(species, s)
    {
        Info<< std::setw(15) << s ;
    }
    Info << endl
        << "========================================================================================"
        << "============================================================================" << endl;


    for (int i = 1; i<10000; i++)
    {
        Info<< std::setw(10) << i;
        forAll(species, s)
        {
            Info<<  std::setw(15) << c.at(s);
        }
        /*Info<< endl;
        Info<< std::setw(10) << i;
        forAll(species, s)
        {
            Info<<  std::setw(15) << chem_.calculateOmega(s, T, c);
        }
        */
        Info<< "\n";

        scalar dt=1e-13;
        //- Omega *dt is the change
        forAll(species,s)
        c[s] += chem_.calculateOmega(s, T, c)*dt; 

        forAll(species, s)
        {
            if (c[s] < 0) c[s] = 0;
        }

    }
    //- Start solving the ODE
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
