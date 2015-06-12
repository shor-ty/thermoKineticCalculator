/*---------------------------------------------------------------------------*\
| =========                |                                                  |
| \\      /  A utomatic    | Holzmann-cfd                                     |
|  \\    /   F lamelet     | Version: 1.0                                     |
|   \\  /    C onstructor  | Web: www.Holzmann-cfd.de                         |
|    \\/                   |                                                  |
\*---------------------------------------------------------------------------*/
/*
»
»
»
»
»
»
»
»
»
»
\*---------------------------------------------------------------------------*/
//- system headers
#include <fstream>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>

//- user def. headers
#include "functions.hpp"
#include "../definitions/IOStream.hpp"
#include "../species/species.hpp"
#include "../database/elements.hpp"
#include "../reactions/reactions.hpp"

//- file functions

    //- open file and return the whole file as string array
    const stringField openFile(const normalString& fileName)
    {
        std::ifstream file;

        //- open file
        file.open(fileName.c_str(), std::ios::in);

        //- check file
        if(!file.good())
        {
            std::cerr<< "\n ++ ERROR: Could not open file \""
                     << fileName << "\"..."
                     << "\n ++ Error occur in file " << __FILE__
                     << " line " << __LINE__ << std::endl;
            std::terminate();
        }

        //- stored line
        normalString fileLine;

        //- file content
        stringField fileContent(0);

        //- read the whole file
        while(!file.eof())
        {
            std::getline(file, fileLine);
            fileContent.push_back(fileLine);
        }
        return fileContent;
    }

    //- read chemical kinetic file
    void readChemKinThermo
    (
        const normalString& fileChemKin,
        const normalString& fileThermo,
        std::vector<Species>& species,
        std::vector<Reactions>& reactions
    )
    {
        stringField fileContent = openFile(fileChemKin);
        int countI{0};

        forAll(fileContent, line)
        {
            //- ELEMENTS
            if (fileContent[line] == "ELEMENTS")
            {
                for(;;line++)
                {
                    //- Do nothing

                    if (fileContent[line] == "END")
                    {
                        line++;
                        break;
                    }
                }
            }

            //- SPECIES
            if (fileContent[line] == "SPECIES")
            {
                for(;;line++)
                {
                    std::istringstream tmp(fileContent[line]);

                    std::vector<std::string> species_
                    {
                        std::istream_iterator<std::string>{tmp},
                        std::istream_iterator<std::string>{}
                    };

                    //- Each species eq. one object
                    forAll(species_, i)
                    {
                        species.resize(species.size()+1);
                        species[countI].setName(species_[i]);
                        countI++;
                    }


                    if (fileContent[line] == "END")
                    {
                        //- Remove first and last entry
                        //  + SPECIES
                        //  + END
                        species.erase(species.begin());
                        species.erase(species.end());

                        //- At least calculate molecular weight
                        forAll(species,i)
                        {
                            species[i].setMW
                            (
                                calcMolecularWeight(species[i].name())
                            );
                        }

                        line++;
                        break;
                    }
                }
            //- SPECIES
            }

            //- reset count variable for reaction count
            countI = 0;

            //- reactions
            if (fileContent[line] == "REACTIONS")
            {
                for(;;line++)
                {
                    //-
                    if
                    (
                        fileContent[line] != "REACTIONS" &&
                        fileContent[line] != "END" &&
                        !fileContent[line].empty()
                    )
                    {
                        //- split line into array
                        std::istringstream tmp(fileContent[line]);

                        stringField lineContent_
                        {
                            std::istream_iterator<std::string>{tmp},
                            std::istream_iterator<std::string>{}
                        };

                        if (lineContent_[0].find('=') != std::string::npos)
                        {
                            //- add new object
                            reactions.resize(reactions.size()+1);

                            //- set the elementar reaction
                            reactions[countI].setElementarReaction
                            (
                                lineContent_[0]
                            );

                            //- check out if only forward is used
                            //  + sign =>
                            if
                            (
                                lineContent_[0].find('>') != std::string::npos
                                &&
                                lineContent_[0].find('<') == std::string::npos
                            )
                            {
                                reactions[countI].set_kf();
                            }
                            //  + sign <=>
                            else if
                            (
                                lineContent_[0].find('>') != std::string::npos
                                &&
                                lineContent_[0].find('<') != std::string::npos
                            )
                            {
                                reactions[countI].set_kf();
                                reactions[countI].set_kb();
                            }
                            //  + sign <=
                            else if
                            (
                                lineContent_[0].find('>') == std::string::npos
                                &&
                                lineContent_[0].find('<') != std::string::npos
                            )
                            {
                                reactions[countI].set_kb();
                            }
                            //  + sign =
                            else
                            {
                                reactions[countI].set_kf();
                                reactions[countI].set_kb();
                            }

                            //- set variable for arrhenius eqn.
                            reactions[countI].setArrheniusCoeffs
                            (
                                stod(lineContent_[1]),
                                stod(lineContent_[2]),
                                stod(lineContent_[3])
                            );

                            //- check if LOW
                            //  + check next line for LOW
                            std::istringstream tmp2(fileContent[line+1]);

                            stringField lineContent2_
                            {
                                std::istream_iterator<std::string>{tmp2},
                                std::istream_iterator<std::string>{}
                            };

                            if (lineContent2_[0] == "LOW/")
                            {
                                //- set LOW
                                reactions[countI].setLOW();

                                //- set variable for arrhenius eqn.
                                reactions[countI].setArrheniusCoeffs
                                (
                                    stod(lineContent2_[1]),
                                    stod(lineContent2_[2]),
                                    stod(lineContent2_[3])
                                );

                                //- check if TROE
                                //  + check next line for TROE (after LOW)
                                std::istringstream tmp3(fileContent[line+2]);

                                stringField lineContent3_
                                {
                                    std::istream_iterator<std::string>{tmp3},
                                    std::istream_iterator<std::string>{}
                                };

                                if (lineContent3_[0] == "TROE/")
                                {
                                    //- set TROE
                                    reactions[countI].setTROE();

                                    //- check if Tss is used
                                    scalar Tss{0};

                                    if (lineContent3_.size() == 6)
                                    {
                                        Tss = stod(lineContent3_[4]);
                                    }

                                    //- set TROE coeffs
                                    reactions[countI].setTROECoeffs
                                    (
                                        stod(lineContent3_[1]),
                                        stod(lineContent3_[2]),
                                        stod(lineContent3_[3]),
                                        Tss
                                    );
                                }
                            }
                            countI++;
                        }
                    }

                    if (fileContent[line] == "END")
                    {
                        line++;
                        break;
                    }
                }
            }
        //- chemKin
        }

        //- Get thermodynamic data from NASA polynomials for each species
        fileContent = openFile(fileThermo);
        bool skip{true};
        int id{-1};

        //- File description from chemkin
        //-------------------------------------------------------------------------------------------------
        //- Line Number         Content                                     Format          Column
        //-------------------------------------------------------------------------------------------------
        //-     1               THERMO (or THERMO ALL »a«)                  Free            Any
        //-     2»b«    Temperatur ranges for 2 sets of coefficents:        Float           1-30
        //-             lowest, common and highest
        //-     3       Species name (must start in column 1)               Char            1-18
        //-             Date (not used                                      Char            19-24
        //-             Atomic symbols and formula                          Char/Int        25-44
        //-             Phase of species (S,L or G)                         Char            45
        //-             Low temperature                                     Float           46-55
        //-             High temperature                                    Float           56-65
        //-             Common temperature                                  Float           66-73
        //-             Atomic symbols and formula (if needed, else blank)  Char/Int        74-78
        //-             The integer 1                                       Int             80
        //-             Atomic symbols and formula (if needed, else blank)  Char/Int        81-100
        //-     4       Coefficnents a1 - a5 (for eqn1-3)                   double          1-75
        //-             for upper temperature interval
        //-             The integer 2                                       Int             80
        //-     5       Coefficnents a6, a7 for upper temperature           double          1-75
        //-             and a1 - a3 for lower temperature interval
        //-             The integer 3                                       Int             80
        //-     6       Coefficnents a4 - a7 for lower temperature          double          1-60
        //-             interval
        //-             The integer 4                                       Int             80
        //-             Repeat lines 3 - 6 for each species
        //-             End (optional, end of thermodynamic data)           Free            Any
        //-------------------------------------------------------------------------------------------------
        //-    »a«      FIXME
        //-    »b«      FIXME
        //- FIXME       Additionall for more accuracy!
        //-------------------------------------------------------------------------------------------------
        //-
        //-     Eqn1:   C_p^0/R     = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
        //-     Eqn2:   H^0/(R*T)   = a1 + a2/2*T^2 + a3/3*T^2 + a4/4*T^3 + a5/5*T^4 + a6/T
        //-     Eqn3:   S^0/R       = a1*ln(T) + a2*T + a3/2*T^2 + a4/3*T^3 + a5/4*T^4 +a7
        //-     Eqn4:   G^0         = H^0 - T*S^0
        //-------------------------------------------------------------------------------------------------

        //- temp polyCoeffs field
        scalarField polyCoeffs;

        forAll(fileContent, line)
        {
            //- analyse first line
            //- FIXME - include "THERMO ALL" extension
            if(line == 0 && (fileContent[line] != "THERMO"))
            {
                std::cerr<< "\n ++ ERROR: first line in thermodynamic file. "
                         << "THERMO\" not found in first line..."
                         << "\n ++ Error occur in file " << __FILE__
                         << " line " << __LINE__ << std::endl;
                std::terminate();
            }

            //- FIXME LINE 2

            //- polyCoeffs after line 2 till end
            if(line > 1)
            {
                //- LINE OF THE INTEGER 1
                if(fileContent[line][79] == '1')
                {
                    //- check if species are used
                    forAll(species, speciesI)
                    {
                        normalString tmp = fileContent[line].substr(0,18);
                        tmp.erase
                        (
                            std::remove(tmp.begin(), tmp.end(), ' '),
                            tmp.end()
                        );

                        if (species[speciesI].name() == tmp)
                        {
                            skip = false;
                            id = speciesI;
                        }
                    }

                    //- get all data from first line
                    if (skip == false)
                    {
                        //- set thermodynamic bool to true
                        species[id].thermodynamicTrue();

                        scalar lowTemp = stod(fileContent[line].substr(45,54));
                        scalar comTemp = stod(fileContent[line].substr(65,72));
                        scalar higTemp = stod(fileContent[line].substr(55,64));

                        species[id].setPolyTemperature
                        (
                            lowTemp,
                            comTemp,
                            higTemp
                        );

                        //- set phase
                        species[id].setPhase
                        (
                            fileContent[line].substr(45,1)
                        );
                    }
                }
                //- LINE OF THE INTEGER 2, 3, 4
                else if
                (
                    (
                        fileContent[line][79] == '2' ||
                        fileContent[line][79] == '3' ||
                        fileContent[line][79] == '4'
                    ) &&
                    skip == false
                )
                {
                    forAll(fileContent[line], pos)
                    {
                        if(pos <= 60)
                        {
                            polyCoeffs.push_back
                            (
                                std::stod(fileContent[line].substr(pos,pos+15))
                            );
                            pos+=14;
                        }
                    }

                    //- INTEGER 4 -> store polyCoeffs
                    if(fileContent[line][79] == '4')
                    {
                        species[id].setPolyCoeffs(polyCoeffs);
                        polyCoeffs.clear();
                    }

                    //- if last line of polynomials reached, set skip to true again
                    if (fileContent[line][79] == '4')
                    {
                        skip = true;
                    }
                }
            }
        //- thermodynamic
        }

        //- check thermodynamic status
        //  + each species need thermodynamics
        //  + if failed - abort
        forAll(species,i)
        {
            if (!species[i].thermodynamicStatus())
            {
                std::cerr<< "\n ++ ERROR: For species "
                         << species[i].name()
                         << " no thermodynamic data found, no NASA polynoms"
                         << std::endl;
                std::terminate();
            }

        }
    }

    //- read afcDict
    void readAFCDict
    (
        const normalString& fileAFCDict,
        std::vector<Species>& species
    )
    {
        stringField fileContent = openFile(fileAFCDict);
        bool fuelFound{false};
        bool oxidizerFound{false};

        forAll(fileContent, line)
        {

            normalString lineContent = fileContent[line];

            //- remove whitespaces
            removeSpace(lineContent);


            if
            (
                lineContent == "moleFractionFuel" ||
                lineContent == "moleFractionOxidizer" ||
                fuelFound == true ||
                oxidizerFound == true
            )
            {
                //- check next line and set bool to true (now in dictionary)
                if (!fuelFound && !oxidizerFound)
                {
                    normalString tmp = fileContent[line+1];
                    removeSpace(tmp);

                    if ( tmp == "{")
                    {
                        if (lineContent == "moleFractionFuel")
                        {
                            fuelFound = true;
                        }
                        if (lineContent == "moleFractionOxidizer")
                        {
                            oxidizerFound = true;
                        }
                    }
                    else
                    {
                        std::cerr<< "\n ++ ERROR: after moleFraction<F/O>"
                                 << " >> { << missing, no dictionary found"
                                 << "\n ++ Error occur in file " << __FILE__
                                 << " line " << __LINE__ << std::endl;
                        std::terminate();
                    }
                    line++;
                }
                //- dictionary end reached "}"
                if (lineContent == "}")
                {
                    fuelFound = false;
                    oxidizerFound = false;
                }

                //- get species of fuel
                if (fuelFound || oxidizerFound)
                {
                    //- re-read due to removed whitespaces
                    //  TODO, change this
                    lineContent = fileContent[line];

                    //- split line into array
                    std::istringstream tmp(lineContent);

                    stringField lineContent_
                    {
                        std::istream_iterator<std::string>{tmp},
                        std::istream_iterator<std::string>{}
                    };

                    //- find species
                    forAll(species, id)
                    {
                        if (species[id].name() == lineContent_[0])
                        {
                            if (fuelFound)
                            {
                                species[id].setFuel();
                            }
                            if (oxidizerFound)
                            {
                                species[id].setOxidizer();
                            }
                            species[id].setMolFraction
                            (
                                stod(lineContent_[1])
                            );
                        }
                    }

                }



            }
        }

        //- check if X = 1
        scalar X{0};
        std::cout<< "Checking chemical composition of fuel:\n";
        forAll(species, id)
        {
            if (species[id].fuel())
            {
                std::cout<< " ++ " << species[id].name() << "\t X = " << species[id].molFraction() << "\n";
                X += species[id].molFraction();
            }
        }
        std::cout<< "Fuel mol fraction = " << X << "\n";

        //- if sum <> 1
        if (X != 1)
        {
            std::cerr<< "\n" << "Fuel mol fraction "
                     << "is not equal to 1\n"
                     << "\n ++ Error occur in file " << __FILE__
                     << " line " << __LINE__ << std::endl;
        }

        //- reset X
        X = 0;
        std::cout<< "Checking chemical composition of oxidizer:\n";
        forAll(species, id)
        {
            if (species[id].oxidizer())
            {
                std::cout<< " ++ " << species[id].name() << "\t X = " << species[id].molFraction() << "\n";
                X =+ species[id].molFraction();
            }
        }
        std::cout<< "Oxidizer mol fraction = " << X << "\n";

        //- if sum <> 1
        if (X != 1)
        {
            std::cerr<< "\n" << "Oxidizer mol fraction "
                     << "is not equal to 1\n"
                     << "\n ++ Error occur in file " << __FILE__
                     << " line " << __LINE__ << std::endl;
        }
    }

    //- remove whitespace from string
    void removeSpace(normalString& str)
    {
        str.erase
        (
            std::remove_if
            (
                str.begin(),
                str.end(),
                [](char c)
                {
                    return (c =='\r' || c =='\t' || c == ' ' || c == '\n');
                }
            ),
            str.end()
        );
    }

    //- calculate molecular weight
    scalar calcMolecularWeight(const normalString& species)
    {
        //- class of molecular weights
        Elements element;

        //- molecular weight [mol/kg]
        scalar mW{0};

        //- there should be a better way to perform that stuff | simple c style
        //  get elements and the number

            normalString tmp = species;
            normalString lP, aP, nP;

            //- inert gas
            if (tmp == "HE" || tmp == "AR")
            {
                mW = element.atomicWeight(tmp);
            }
            else
            {
                for (unsigned int i=0; i<tmp.length(); i++)
                {
                    //- actual position
                    aP = tmp[i];

                    //- next position
                    if (i+1 <= tmp.length())
                    {
                        nP = tmp[i+1];
                    }
                    else
                    {
                        nP = aP;
                    }

                    //- last position
                    if (i > 0)
                    {
                        lP = tmp[i-1];
                    }
                    else
                    {
                        lP = aP;
                    }

                    //- molecule
                    if (aP.find_first_of("0123456789") != std::string::npos)
                    {
                        if
                        (
                            nP.find_first_of("0123456789") != std::string::npos)
                        {
                            mW =+ element.atomicWeight(lP)
                                * stod(tmp.substr(i,i+1));

                            i++;
                        }
                        else
                        {
                            mW =+ element.atomicWeight(lP)
                                * stod(aP);
                        }
                    }
                    //- element
                    else
                    {
                        mW = element.atomicWeight(aP);
                    }
                }
            }

            //- return the molecular weight of species
            return mW;
    }


    //- calculate stoichiometric mixture fraction Zst
    scalar stochiometricMF(const std::vector<Species>& species)
    {
        scalar Zst{0};

        return Zst;
    }
