/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2015 Holzmann-cfd
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Creator.

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
#include "interpreter.hpp"
#include "scalar.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "chemistry.hpp"
#include "thermo.hpp"
#include "transport.hpp"
#include "properties.hpp"
//#include "mixtureFraction.hpp"
//#include "numerics.hpp"
//#include "LUDecompose.hpp"

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

    Info<< Header() << endl;

    string file_AFC;
    string file_Thermo;
    string file_Transport;
    string file_Chemistry;

    //- Build Interpreter object which makes a summary at the end if the switch
    //  is set on (-interprete)
    Interpreter interpreter;

    //- Check if all arguments are correct 
    if (argc != 3)
    {
        ErrorMsg
        (
            "    Program needs two arguments.\n\n"
            "    For calculating flamelets you have to use\n"
            "      ./automaticFlameletCreator -AFCDict $pathToFile\n\n",
            __FILE__,
            __LINE__
        ); 
    }
    else
    {
        file_AFC = string(argv[2]);
    }

    //- Temporary Reader to get the file paths
    {
        PropertiesReader tmp (file_AFC);
        file_Thermo = tmp.path("thermodynamic");
        file_Transport = tmp.path("transport");
        file_Chemistry = tmp.path("chemistry");
    }

    Thermo thermo(file_Thermo);

    Transport transport(file_Transport, thermo);

    Chemistry chemistry(file_Chemistry, thermo);

    Properties properties(file_AFC, thermo, chemistry);

    Info<< " c-o All data read successfully\n" << endl;

    Info<< " c-o Check transportData (all species available)\n" << endl;

    //- Proof if all species used in chemistryData have transportData
    {
       const wordList& speciesCH = chemistry.species();
       const wordList& speciesTD = transport.species();

       forAll(speciesCH, i)
       {
           bool found{false};

           forAll(speciesTD, j)
           {
               if (i == j)
               {
                   found = true;
               }
           }

           //- Chemistry species not found in transportData
           if (!found)
           {
               ErrorMsg
               (
                   "Species " + i + " was not found in the transportData class."
                   "\n    Please be sure that the transport properties of "
                   "this species is available in the transport file.",
                   __FILE__,
                   __LINE__
               );
           }
       }

        //- Insert chemical species to transportData
        transport.insertChemistrySpecies(speciesCH);
    }

    Info<< " c-o Check thermodynamicData (all species available)\n" << endl;

    //- Proof if all species used in chemistryData have thermodynamicData
    {
        const wordList& speciesCH = chemistry.species();
        const wordList& speciesTHD = thermo.species();

        forAll(speciesCH, i)
        {
            bool found{false};

            forAll(speciesTHD, j)
            {
                if (i == j)
                {
                    found = true;
                }
            }

            //- Chemistry species not found in thermodynamicData
            if (!found)
            {
                ErrorMsg
                (
                    "Species " + i + " was not found in the thermodynamicData"
                    "class.\n    Please be sure that the thermodynamic propert"
                    "ies of this species is available in the termo file.",
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
            speciesAFC.push_back(i);
        }

        //- Check if species are in chemistry
        const wordList& speciesCH = chemistry.species();

        forAll(speciesAFC, i)
        {
            bool found{false};

            forAll(speciesCH, j)
            {
                if (i == j)
                {
                    found = true;
                }
            }

            //- afc species not found in chemistry, thermo or transport
            if (!found)
            {
                ErrorMsg
                (
                    "Species " + i + " was not found in the chemistryData "
                    "class.\n    Please be sure that the species is available "
                    "in all files.",
                    __FILE__,
                    __LINE__
                );
            }
        }

        Info<< " c-o Data O.K.\n" << endl;
    }
    
    //- Interprete data 
    if (interpreter.analyze())
    {
        Info<< " c-o Interprete data ...\n" << endl;

        transport.fitCurves();

        interpreter.summary(chemistry, thermo, transport);

        Footer(startTime);

        return 0;
    }

    //- Fit data to polynomials
    transport.fitCurves();


//- TH::200217 
//
//    //- Now we can proceed doing some other stuff after the check
//    {
//        //- Calculate adiabatic enthalpy of fuel and oxidizer
//        properties.calcProperties();
//
//        //- Calculate adiabatic flame temperature (simple estimate)
//        properties.calcAdiabaticTemperature();
//
//        transport.calcAtomComposition();
//    }
//
//    //- Calculate first flamelet (initial - equilibrium)
//    AdiabaticFlamelet adiabaticFlamelet;
//
//    //- First lookUpTable for first scalar dissipation rate
//    //  and adiabatic enthalpy (no defect)
//    /*{
//        Info<< "    ... pre-calculation (get it burn)\n" << endl;
//
//        //- First scalar dissipation rate TODO have to be sorted before
//        const scalar& sDR1 = properties.sDRs(0);
//
//        // - Discretisation of mixture fraction Z
//        const unsigned int& Zpoints = properties.nZPoints();
//
//        const scalar delta = 1./Zpoints;
//
//        Info<< "    ... create adiabatic flamelet for " << sDR1
//            << " Hz\n" << endl;
//
//        //- All points in adiabatic flamelet
//        for(unsigned int i=0; i <= Zpoints; i++)
//        {
//            scalar zPointValue = i*delta;
//
//            adiabaticFlamelet.push_back
//            (
//                MixtureFraction
//                (
//                    chemistry,
//                    properties,
//                    thermo,
//                    transport,
//                    zPointValue,
//                    scalar(0)
//                )
//            );
//        }
//
//        Footer(startTime);
//
//        return 0;
//    }*/
//
//    Info<< " c-o Overview of Look-Up-Tables ...\n" << endl;
//
//    //- Definition
//    //  |
//    //  |-> scalarDissipationRate
//    //      |
//    //      |-> MixtureFraction point [0-1]
//    //          |
//    //          |-> Class of MixtureFraction
//    //map<word, map<unsigned int, MixtureFraction> > flame;
//
//    lookUpTable lookUpTables;
//
//    //- The first flamelet is always complete in equilibrium (sDR -> very low)
//    //  and has no enthalpy defect (adiabatic flamelet)
//    MixtureFraction adiabaticFlamelet
//    (
//        chemistry,
//        thermo,
//        transport,
//        properties,
//        scalar(1e-6),
//        scalar(0)
//    );
//
//    /*{
//        const scalarField defects = properties.defects();
//        const scalarField sDRs = properties.sDRs();
//
//        //- Adiabatic flamelet (always available)
//        {
//            Info<< "    ... for adiabatic condition\n\n";
//
//            lookUpTables.push_back(vector<MixtureFraction>());
//
//            forAll(sDRs, rate)
//            {
//                Info<< "        Create flamelets for scalar dissipation rate "
//                    << rate << " Hz\n";
//
//                lookUpTables[0].push_back
//                (
//                    MixtureFraction
//                    (
//                        chemistry,
//                        thermo,
//                        transport,
//                        properties,
//                        rate,
//                        scalar(0)
//                    )
//                );
//
//            //- sDR 
//            }
//            Info<< "\n";
//
//        //- adiabatic condition
//        }
//    
//        //- Counter (0 is for adiabatic)
//        unsigned int c{1};
//        
//        //- Non-adiabatic flamelets
//        forAll(defects, defect)
//        {
//            if (defect == 0)
//            {
//                break;
//            }
//
//            Info<< "    ... for enthalpy defect " << defect << " J/kg\n\n";
//
//            lookUpTables.push_back(vector<MixtureFraction>());
//
//            forAll(sDRs, rate)
//            {
//                Info<< "        Create flamelets for scalar dissipation rate "
//                    << rate << " Hz\n";
//
//                lookUpTables[c].push_back
//                (
//                    MixtureFraction
//                    (
//                        chemistry,
//                        thermo,
//                        transport,
//                        properties,
//                        rate,
//                        defect
//                    )
//                );
//
//            //- sDR 
//            }
//            Info<< "\n";
//
//        //- defect 
//        }
//        Info<< "\n";
//    }*/
//
//    //- Calculate initial solution
//    Numerics num(chemistry);
//
//    Info<< " c-o Calculate initial solution for adiabatic flamelet\n\n";
//    {
//        num.solveForInitialSolution(adiabaticFlamelet);
//
//        Info<< "\n";
//    }
//
//    Info<< " c-o Get it burn - calculate adiabatic flamelet\n\n";
//    {
//        num.solveAdiabaticFlamelet(adiabaticFlamelet);
//    }
//
//    //- Calculation start
//    /*Info<< " c-o Start flamelet calculation\n" << endl;
//    {
//        const scalarField sDRs = properties.sDRs();
//        const scalarField defects = properties.defects();
//        Numerics num;
//
//        //- Defect loop
//        forAll(defects, defect)
//        {
//            Info<< "    c-o Calculate Look-Up-Table with defect: " << defect
//                << " J/kg\n\n";
//
//            //- Scalar dissiaption loop
//            forAll(sDRs, rate)
//            {
//                Info<< "      c-o Calculate flamelet for " << rate << " Hz\n"
//                    << "      --------------------------------------------\n";
//
//                MixtureFraction& flamelet = lookUpTables[defect][rate];
//
//                scalar saveTime{0};
//                
//                //- Time loop
//                for
//                (
//                    scalar runTime = 0;
//                    runTime < properties.runTime();
//                )
//                {
//                    //- Write first 
//                    if (runTime == 0)
//                    {
//                        flamelet.write();
//                    }
//                    //- Runtime
//                    runTime += properties.deltat();
//
//                    //- Update time
//                    properties.updateCurrentTime(runTime);
//
//                    //- (SaveTime) increase
//                    saveTime += properties.deltat();;
//
//                    Info<< "        Time: " << runTime << " s\n";
//
//                    num.solveFlamelet(flamelet, rate, properties.deltat());
//
//                    //- Simple save algorithm
//                    if (saveTime >= properties.write())
//                    {
//                        flamelet.write();
//
//                        saveTime = 0;
//                    }
//                }
//            }
//        }
//    }*/
//
    Footer(startTime);

    return 0;
}


// ************************************************************************* //

