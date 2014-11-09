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
#include <iostream>

//- user def. headers
#include "functions.hpp"


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
            std::cerr << "\n ++ ERROR: Could not open file \"" << file << "\"...";
            std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__ << std::endl;
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

    //- create the thermodynamic objects and store all files
    std::vector<Thermodynamic> createThermodynamicObjects(const stringField& thermodynamicFileContent)
    {
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

        //- dynamic thermodynamic vector class
        std::vector<Thermodynamic> thermo;

        //- temp polyCoeffs
        scalarField polyCoeffs_;

        forAll(thermodynamicFileContent, line)
        {
            //- analyse first line
            //- FIXME - include "THERMO ALL" extension
            if(line == 0 && (thermodynamicFileContent[line] != "THERMO"))
            {
                std::cerr << "\n ++ ERROR: first line in thermodynamic file. \"THERMO\" not found in first line...";
                std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__ << std::endl;
                std::terminate();
            }

            //- FIXME LINE 2

            //- polyCoeffs after line 2 till end
            if(line > 1)
            {
            //- LINE OF THE INTEGER 1
                if(thermodynamicFileContent[line][79] == '1')
                {
                    //- new entry
                    thermo.resize(thermo.size()+1);

                    //- temp temperature bounds
                    scalarField temperatureBound_(3);
                    temperatureBound_[0] = std::stod(thermodynamicFileContent[line].substr(46,10));
                    temperatureBound_[2] = std::stod(thermodynamicFileContent[line].substr(56,10));
                    temperatureBound_[1] = std::stod(thermodynamicFileContent[line].substr(66,10));

                    //- set temperature bounds
                    thermo[thermo.size()-1].setTemperatureBounds(temperatureBound_);

                }
                //- LINE OF THE INTEGER 2, 3, 4
                else if((thermodynamicFileContent[line][79] == '2' || thermodynamicFileContent[line][79] == '3' || thermodynamicFileContent[line][79] == '4'))
                {
                    forAll(thermodynamicFileContent[line], pos)
                    {
                        if(pos <= 60)
                        {
                            polyCoeffs_.push_back(std::stod(thermodynamicFileContent[line].substr(pos,pos+15)));
                            pos+=14;
                        }
                    }

                    //- INTEGER 4 -> store polyCoeffs
                    if(thermodynamicFileContent[line][79] == '4')
                    {
                        thermo[thermo.size()-1].setPolyCoeffs(polyCoeffs_);
                        polyCoeffs_.clear();
                    }
                }
            }


            //- FIXME extend to
            //if()
        }            //forAll(polyCoeffs_, i) std::cout << polyCoeffs_[i] << "\n";
        return thermo;
    }

