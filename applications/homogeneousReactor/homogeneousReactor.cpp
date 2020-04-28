/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     T hermo       | Open Source Thermo-Kinetic Library
  c-o-o-c     K iknetic     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Creator.

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

Description

    Transient thermo-kinetic 0D calculator for detailed chemistry analysis.


\*---------------------------------------------------------------------------*/

#include "definitions.hpp"
#include "idealReactorProperties.hpp"
#include "transport.hpp"
#include "thermo.hpp"
#include "chemistry.hpp"
#include "interpreter.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace TKC;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char** argv)
{

    const std::clock_t startTime = clock();

    Info<< Header() << endl;

    //- Create Objects for calculation
    IdealReactorProperties properties("");

    Thermo thermo(properties.thermo());

    Transport transport(properties.transport(), thermo);

    Chemistry chemistry(properties.chemistry(), thermo);

    //- Interprete data and store for analysis in files
    if (properties.interprete())
    {
        Interpreter interpreter;

        interpreter.summary(transport, thermo, chemistry);

        Footer(startTime);
        return 0;
    }



    //- Initialize 0D calculation


    //- Solve chemistry and save data regarding the use input



    Footer(startTime);

    return 0;
}


// ************************************************************************* //

