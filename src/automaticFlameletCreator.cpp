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

Description


\*---------------------------------------------------------------------------*/

#include "typedef.hpp" 
#include "chemistry.hpp"
#include "thermo.hpp"
#include "transport.hpp"
#include "mixtureFraction.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace AFC;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main
(
    int argc,
    char** argv
)
{
    Header();

    string file_AFC;
    string file_Thermo;
    string file_Transport;
    string file_Chemistry;

    stringList test;

    //- Arguments
    if (argc != 9)
    {
        FatalError
        (
            "    Program needs eight arguments.\n\n"
            "    ./automaticFlameletCreator\n"
            "      -transport $pathToFile\n"
            "      -thermodynamic $pathToFile\n"
            "      -chemistry $pathToFile\n"
            "      -AFCDict $pathToFile\n",
            __FILE__,
            __LINE__
        ); 
    }
    else
    {
        for (int i=0; i<argc; i++)
        {
            if (string(argv[i]) == "-transport")
            {
                file_Transport = string(argv[i+1]);
            }
            else if (string(argv[i]) == "-thermodynamic")
            {
                file_Thermo = string(argv[i+1]);
            }
            else if (string(argv[i]) == "-chemistry")
            {
                file_Chemistry = string(argv[i+1]);
            }
            else if (string(argv[i]) == "-AFCDict")
            {
                file_AFC = string(argv[i+1]);
            }
        }

        if
        (
            file_Transport.empty()
         || file_Thermo.empty()
         || file_Chemistry.empty()
         || file_AFC.empty()
        )
        {
            FatalError
            (
                "    Tranport, Thermo, Chemistry or AFCDict is missing.\n"
                "    Check the command that you use for running the "
                "application.",
                __FILE__,
                __LINE__
            );
        }
    }

    Chemistry chemistry(file_Chemistry);

    Thermo thermo(file_Thermo, chemistry.thermo());

    Transport transport(file_Transport);

    MixtureFraction Z0(file_AFC);

    Info<< " c-o All data read successfully\n" << endl;

    Info<< " c-o Re-organize data structure\n" << endl;

    return 0;
}


// ************************************************************************* //
