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
#include "properties.hpp"

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

    Properties properties(file_AFC);

    Info<< " c-o All data read successfully\n" << endl;

    Info<< " c-o Check transportData (all species available)\n" << endl;

    //- Proof if all species used in chemistryData have transportData
    {
        wordList speciesCH = chemistry.species();
        wordList speciesTD = transport.species();

        forAll(speciesCH, i)
        {
            bool found{false};

            forAll(speciesTD, j)
            {
                if (speciesCH[i] == speciesTD[j])
                {
                    found = true;
                }
            }

            //- Chemistry species not found in transportData
            if (!found)
            {
                FatalError
                (
                    "    Species " + speciesCH[i] + " was not found in "
                    "the transportData class.\n"
                    "    Please be sure that the transport properties of "
                    "this species is available\n"
                    "    in the transport file.",
                    __FILE__,
                    __LINE__
                );
            }
        }
    }

    Info<< " c-o Check thermodynamicData (all species available)\n" << endl;

    //- Proof if all species used in chemistryData have thermodynamicData
    {
        wordList speciesCH = chemistry.species();
        wordList speciesTHD = thermo.species();

        forAll(speciesCH, i)
        {
            bool found{false};

            forAll(speciesTHD, j)
            {
                if (speciesCH[i] == speciesTHD[j])
                {
                    found = true;
                }
            }

            //- Chemistry species not found in thermodynamicData
            if (!found)
            {
                FatalError
                (
                    "    Species " + speciesCH[i] + " was not found in "
                    "the thermodynamicData class.\n"
                    "    Please be sure that the thermodynamic properties of "
                    "this species is available\n"
                    "    in the thermodynamic file.",
                    __FILE__,
                    __LINE__
                );
            }
        }
    }

    Info<< " c-o Check if species used for combustion are available\n" << endl;

    //- Proof if all species in afcDict are in chemistryData
    {
        wordList speciesOAFC = properties.speciesOxidizer();
        wordList speciesFAFC = properties.speciesFuel();
        wordList speciesAFC = speciesOAFC;

        //- Add fuel species
        forAll(speciesFAFC, i)
        {
            speciesAFC.push_back(speciesFAFC[i]);
        }

        wordList speciesCH = chemistry.species();

        forAll(speciesAFC, i)
        {
            bool found{false};

            forAll(speciesCH, j)
            {
                if (speciesAFC[i] == speciesCH[j])
                {
                    found = true;
                }
            }

            //- afc species not found in chemistryData 
            if (!found)
            {
                FatalError
                (
                    "    Species " + speciesAFC[i] + " was not found in "
                    "the chemistryData class.\n"
                    "    Please be sure that the species is available in all "
                    "files; chemistry, tranpsort\n"
                    "    and the thermodynamic file.",
                    __FILE__,
                    __LINE__
                );
            }
        }

        Info<< " c-o Data O.K.\n" << endl;
    }

    return 0;
}


// ************************************************************************* //
