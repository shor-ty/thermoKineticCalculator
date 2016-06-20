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

#include "time.h"
#include "typedef.hpp" 
#include "chemistry.hpp"
#include "thermo.hpp"
#include "transport.hpp"
#include "properties.hpp"
#include "mixtureFraction.hpp"
#include "numerics.cpp"
#include "interprete.cpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace AFC;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main
(
    int argc,
    char** argv
)
{
    const std::clock_t startTime = clock();

    Header();


    string file_AFC;
    string file_Thermo;
    string file_Transport;
    string file_Chemistry;

    //- Arguments
    bool interprete{false};

    if (argc != 9 && argc != 10)
    {
        FatalError
        (
            "    Program needs one or eight arguments.\n\n"
            "    For calculating flamelets you have to use\n"
            "    ./automaticFlameletCreator\n"
            "      -transport $pathToFile\n"
            "      -thermodynamic $pathToFile\n"
            "      -chemistry $pathToFile\n"
            "      -AFCDict $pathToFile\n\n"
            "    For interpreting the data you have to add\n"
            "      -interprete\n",
            __FILE__,
            __LINE__
        ); 
    }
    else
    {
        if (argc == 10)
        {
            interprete = true;
        }

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

    Info<< " c-o Create Objects\n" << endl;

    Chemistry chemistry(file_Chemistry);

    Thermo thermo(file_Thermo, chemistry.thermo());

    Transport transport(file_Transport);

    Properties properties(file_AFC, thermo, chemistry);

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

        //- Insert chemical species to transportData
        transport.chemSpecies(speciesCH);
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

    Info<< " c-o Check if species used in afcDict are available\n" << endl;

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

    //- Now we can proceed doing some other stuff after the check
    {
        //- Calculate adiabatic enthalpy of fuel and oxidizer
        properties.calcProperties();

        //- Calculate adiabatic flame temperature (simple estimate)
        properties.calcAdiabaticTemperature();

        //transport.calcAtomComposition();
    }

    //- Calculate viscosity and binary diffusion coefficients for fitting
    //  procedure
    {
//        transport.prepareFitting(thermo);
    }
    //- Interprete data
    if (interprete)
    {
        Info<< " c-o Interprete data ...\n" << endl;

        AFC::interprete(thermo, chemistry, properties);

        Info<< "\n c-o Interperted all data\n" << endl;

        Footer(startTime);

        return 0;
    }

    //- Calculate first flamelet (initial - equilibrium)
    AdiabaticFlamelet adiabaticFlamelet;

    //- First lookUpTable for first scalar dissipation rate
    //  and adiabatic enthalpy (no defect)
    {
        Info<< "    ... pre-calculation (get it burn)\n" << endl;

        //- First scalar dissipation rate TODO have to be sorted before
        const scalar& sDR1 = properties.sDRs(0);

        // - Discretisation of mixture fraction Z
        const unsigned int& Zpoints = properties.mfPoints();

        const scalar delta = 1./Zpoints;

        Info<< "    ... create adiabatic flamelet for " << sDR1
            << " Hz\n" << endl;

        //- All points in adiabatic flamelet
        for(unsigned int i=0; i <= Zpoints; i++)
        {
            scalar zPointValue = i*delta;

            adiabaticFlamelet.push_back
            (
                MixtureFraction
                (
                    chemistry,
                    properties,
                    thermo,
                    transport,
                    zPointValue,
                    scalar(0)
                )
            );
        }

        Footer(startTime);

        return 0;
    }

    Info<< " c-o Overview of Look-Up-Tables ...\n" << endl;

    //- Definition
    //  |
    //  |-> scalarDissipationRate
    //      |
    //      |-> MixtureFraction point [0-1]
    //          |
    //          |-> Class of MixtureFraction
    //map<word, map<unsigned int, MixtureFraction> > flame;
    LookUpTable lookUpTables;

    {
        const scalarField defects = properties.defects();
        const scalarField sDRs = properties.sDRs();
        const unsigned int Zpoints = properties.mfPoints();

        scalar delta = 1./Zpoints;

        //- Adiabatic flamelet (always available)
        {
            Info<< "    ... for adiabatic condition\n\n";

            lookUpTables.push_back(vector<vector<MixtureFraction> >());

            forAll(sDRs, rate)
            {
                Info<< "        Create flamelets for scalar dissipation rate "
                    << sDRs[rate] << " Hz\n";

                lookUpTables[0].push_back(vector<MixtureFraction>());

                for(unsigned int i=0; i <= Zpoints; i++)
                {
                    scalar zPointValue = i*delta;

                    lookUpTables[0][rate].push_back
                    (
                        MixtureFraction
                        (
                            chemistry,
                            properties,
                            thermo,
                            transport,
                            zPointValue,
                            scalar(0)
                        )
                    );

                //- Z points
                }

            //- sDR points
            }
            Info<< "\n";

        //- adiabatic condition
        }
    
        //- Non-adiabatic flamelets
        forAll(defects, defect)
        {
            Info<< "    ... for enthalpy defect " << defects[defect] << " J/kg\n\n";

            lookUpTables.push_back(vector<vector<MixtureFraction> >());

            forAll(sDRs, rate)
            {
                Info<< "        Create flamelets for scalar dissipation rate "
                    << sDRs[rate] << " Hz\n";

                lookUpTables[defect+1].push_back(vector<MixtureFraction>());

                for(unsigned int i=0; i < Zpoints+1; i++)
                {
                    scalar zPointValue = i*delta;

                    lookUpTables[defect+1][rate].push_back
                    (
                        MixtureFraction
                        (
                            chemistry,
                            properties,
                            thermo,
                            transport,
                            zPointValue,
                            defects[defect]
                        )
                    );

                //- Z points
                }

            //- sDR points
            }
            Info<< "\n";

        //- defect points
        }
        Info<< "\n";
    }

    //- Calculation start
    Info<< " c-o Start flamelet calculation\n" << endl;
    {
        const scalarField sDRs = properties.sDRs();
        const scalarField defects = properties.defects();
        const unsigned int Zpoints = properties.mfPoints();
    
        //- Defect loop
        for
        (
            unsigned int defectNo = 0;
            defectNo < properties.nDefects();
            defectNo++
        )
        {
            Info<< "    Calculate Look-Up-Table with defect: "
                << defects[defectNo]
                << " J/kg\n";

            //- Scalar dissiaption loop
            forAll(sDRs, rate)
            {

                //- Time loo
                for
                (
                    scalar time = 0;
                    time < properties.runTime();
                    time += properties.deltaT()
                )
                {
                    Info<< "Time: " << time << " s\n";
                    //- Equations (in *.pdf file)
                    //  Laminar flamelet model for temperature Eqn (1)
                    //  Laminar flamelet model for species Eqn (2) 
                    calculate
                    (
                        lookUpTables,
                        sDRs[rate],
                        defectNo,
                        Zpoints,
                        properties.deltaT()
                    );
                }
            }
        }
    }

    Footer(startTime);

    return 0;
}


// ************************************************************************* //
