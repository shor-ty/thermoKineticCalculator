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

Class
   AFC::Constants
   
Description
    This class contains all constant data used in the program.

SourceFiles
    constants.cpp

\*---------------------------------------------------------------------------*/

#ifndef Constants_hpp
#define Constants_hpp

#include "typedef.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

namespace Constants
{

/*--------------------------------------------------------------------------*\
                                    Constants
\*--------------------------------------------------------------------------*/

//- Constants and their weight [g/mol]
const map<word, scalar> AW
{
    { "H", 1.00790 },
    { "O", 15.9994 },
    { "C", 12.0110 },
    { "CL", 35.453 },
    { "N", 14.0067 },
    { "AR", 39.9480 },
    { "HE", 4.00260 }
};

//- Converstion from cal to joule [J/cal]
constexpr scalar calToJoule = 4.1858;

//- Converstion from cal to joule [cal/J]
constexpr scalar jouleToCal = 1/calToJoule;

//- Universal gas constant [J/K/mol]
constexpr scalar R = 8.314459848;

//- Universal gas constant [cal/K/mol]
constexpr scalar Rcal = R * jouleToCal;

//- Universal gas constant [erg/K/mol]
constexpr scalar Rerg = 8.314459848e7;

//- Boltzmann Constant [J/K] = [kg m / s^2 K]
constexpr scalar kB = 1.3806485279e-23;

//- Stefan-Boltzmann Constant
//constexpr scalar kSB = 1;

//- Avogadro Constant [1/mol] => R / kB
constexpr scalar N_A = R /kB;

//- Reference pressure (normally 1atm, we use Pa)
constexpr scalar p0 = 1e5;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // Namespace Constants

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // Namespace AFC 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Constants_hpp included

// ************************************************************************* //
