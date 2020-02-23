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

#include "propertiesReader.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::PropertiesReader::PropertiesReader(const string file)
:
    file_(file)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::PropertiesReader::~PropertiesReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AFC::PropertiesReader::read(Properties& data)
{
    Info<< " c-o Reading AFCDict\n" << endl;

    const auto fileContent = readFile(file_);

    //- Reading afcDict
    for (unsigned int line=0; line < fileContent.size(); line++)
    {
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            if (tmp[0] == "pressure")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No pressure specified (" + file_ +")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.insertPressure(stod(tmp[1]));
            }
            else if (tmp[0] == "inertGas")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No inert gas specified (" + file_ +")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.inertSpecies(word(tmp[1]));
            }
            else if (tmp[0] == "fuel")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No fuel specified (" + file_ +")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.fuelSpecies(word(tmp[1]));
            }
            else if (tmp[0] == "oxidizer")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No oxidizer specified (" + file_ +")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.oxidizerSpecies(word(tmp[1]));
            }
            else if (tmp[0] == "mixtureFractionPoints")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No mixtureFractionPoints specified (" + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.insertMFPoints(stoi(tmp[1]));
            }
            else if (tmp[0] == "varianzOfMixtureFractionPoints")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No varianzOfMixtureFractionPoints specified ("
                        + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.insertVMFPoints(stoi(tmp[1]));
            }
            else if (tmp[0] == "enthalpyDefects")
            {
                enthalpyDefects(fileContent, line, data);
            }
            else if (tmp[0] == "scalarDissipationRates")
            {
                scalarDissipationRates(fileContent, line, data);
            }
            else if (tmp[0] == "temperatureOxidizer")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No temperature for oxidizer specified ("
                        + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.insertTemperatureOxidizer(stod(tmp[1]));
            }
            else if (tmp[0] == "temperatureFuel")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No temperature for fuel specified ("
                        + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.insertTemperatureFuel(stod(tmp[1]));
            }
            else if (tmp[0] == "molFractionOxidizer")
            {
                molFractionOxidizer(fileContent, line, data);

                //- Input is mol fraction
                input("mol", data);
            }
            else if (tmp[0] == "molFractionFuel")
            {
                molFractionFuel(fileContent, line, data);
            }
            else if (tmp[0] == "massFractionOxidizer")
            {
                massFractionOxidizer(fileContent, line, data);

                //- Input is mass fraction
                input("mass", data);
            }
            else if (tmp[0] == "massFractionFuel")
            {
                massFractionFuel(fileContent, line, data);
            }
            else if (tmp[0] == "afcControl")
            {
                control(fileContent, line, data);
            }
            else if (tmp[0] == "interpreter")
            {
                interpreter(fileContent, line, data);
            }
        }
    }

    data.check();
}


AFC::string AFC::PropertiesReader::path(const string key)
{
    const auto fileContent = readFile(file_);

    //- Reading afcDict
    for (unsigned int line=0; line < fileContent.size(); line++)
    {
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            if (tmp[0] == key)
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No " + key + " path specified (" + file_ +")",
                        __FILE__,
                        __LINE__
                    );
                }

                return tmp[1];
            }
        }
    }

    //- Suppressing compiler warning warning
    return "";
}


bool AFC::PropertiesReader::interprete()
{
    const auto fileContent = readFile(file_);

    //- Reading afcDict
    for (unsigned int line=0; line < fileContent.size(); line++)
    {
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            if (tmp[0] == "interprete")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "Keyword 'interprete' not specified (" + file_ +")",
                        __FILE__,
                        __LINE__
                    );
                }

                if (tmp[1] == "true" || tmp[1] == "yes") { return true; }
            }
        }
    }

    return false;
}

// * * * * * * * * * * * * * * Helper functions  * * * * * * * * * * * * * * //

void AFC::PropertiesReader::findKeyword
(
    int& start,
    unsigned int& end,
    const stringList& fileContent,
    unsigned int& line
)
{
    for (; line < fileContent.size(); line++)
    {
        //- Split string; delimiter ' ' 
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            if (tmp[0] == "{")
            {
                start = line; 
            }

            //- Search end after keyword found
            if (tmp[0] == "}" && start != -1)
            {
                end = line;
            }

            //- Both found exit loop
            if
            (
                start != 0
             && end != 0
            )
            {
                break;
            }
        }
    }
}


// * * * * * * * * * * * * Data manipulation functions * * * * * * * * * * * //

void AFC::PropertiesReader::enthalpyDefects
(
    const stringList& fileContent,
    unsigned int line,
    Properties& data
)
{
    int dictBegin{-1};
    unsigned int dictEnd{0};
    
    findKeyword(dictBegin, dictEnd, fileContent, line);

    line = dictBegin+1;

    for(; line < dictEnd; line++)
    {
        //- Split string; delimiter ' ' 
        stringList tmp = splitStrAtWS(fileContent[line]);
        
        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            //- If more elements in one line
            if (tmp.size() > 1)
            {
                forEach(tmp, element)
                {
                    data.insertEnthalpyDefects(stod(tmp[element]));
                }
            }
            else
            {
                data.insertEnthalpyDefects(stod(tmp[0]));
            }
        }
    }
}


void AFC::PropertiesReader::scalarDissipationRates
(
    const stringList& fileContent,
    unsigned int line,
    Properties& data
)
{
    int dictBegin{-1};
    unsigned int dictEnd{0};
    
    findKeyword(dictBegin, dictEnd, fileContent, line);

    line = dictBegin+1;

    for(; line < dictEnd; line++)
    {
        //- Split string; delimiter ' ' 
        stringList tmp = splitStrAtWS(fileContent[line]);
        
        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            //- If more elements in one line
            if (tmp.size() > 1)
            {
                forEach(tmp, element)
                {
                    data.insertScalarDissipationRates(stod(tmp[element]));
                }
            }
            else
            {
                data.insertScalarDissipationRates(stod(tmp[0]));
            }
        }
    }
}


void AFC::PropertiesReader::molFractionOxidizer
(
    const stringList& fileContent,
    unsigned int line,
    Properties& data
)
{
    int dictBegin{-1};
    unsigned int dictEnd{0};
    
    findKeyword(dictBegin, dictEnd, fileContent, line);

    line = dictBegin+1;

    for(; line < dictEnd; line++)
    {
        //- Split string; delimiter ' ' 
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            //- If more elements in one line
            if (tmp.size() > 2)
            {
                forEach(tmp, element)
                {
                    unsigned int modul = element%2;

                    if (modul)
                    {
                        if (tmp[element+1].empty())
                        {
                            ErrorMsg
                            (
                                "Problem in dictionary molFractionOxidizer "
                                "(" + file_ + ")",
                                __FILE__,
                                __LINE__
                            );
                        }

                        data.insertCompositionOxidizerMol
                        (
                            tmp[element],
                            stod(tmp[element+1])
                        );
                    }
                }
            }
            else
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "Problem in dictionary molFractionOxidizer "
                        "(" + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.insertCompositionOxidizerMol
                (
                    tmp[0],
                    stod(tmp[1])
                );
            }
        }
    }
}


void AFC::PropertiesReader::massFractionOxidizer
(
    const stringList& fileContent,
    unsigned int line,
    Properties& data
)
{
    int dictBegin{-1};
    unsigned int dictEnd{0};
    
    findKeyword(dictBegin, dictEnd, fileContent, line);

    line = dictBegin+1;

    for(; line < dictEnd; line++)
    {
        //- Split string; delimiter ' ' 
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            //- If more elements in one line
            if (tmp.size() > 2)
            {
                forEach(tmp, element)
                {
                    unsigned int modul = element%2;

                    if (modul)
                    {
                        if (tmp[element+1].empty())
                        {
                            ErrorMsg
                            (
                                "Problem in dictionary massFractionOxidizer "
                                "(" + file_ + ")",
                                __FILE__,
                                __LINE__
                            );
                        }

                        data.insertCompositionOxidizerMass
                        (
                            tmp[element],
                            stod(tmp[element+1])
                        );
                    }
                }
            }
            else
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "Problem in dictionary massFractionOxidizer "
                        "(" + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.insertCompositionOxidizerMass
                (
                    tmp[0],
                    stod(tmp[1])
                );
            }
        }
    }
}


void AFC::PropertiesReader::molFractionFuel
(
    const stringList& fileContent,
    unsigned int line,
    Properties& data
)
{
    int dictBegin{-1};
    unsigned int dictEnd{0};
    
    findKeyword(dictBegin, dictEnd, fileContent, line);

    line = dictBegin+1;

    for(; line < dictEnd; line++)
    {
        //- Split string; delimiter ' ' 
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            //- If more elements in one line
            if (tmp.size() > 2)
            {
                forEach(tmp, element)
                {
                    unsigned int modul = element%2;

                    if (tmp[element+1].empty())
                    {
                        ErrorMsg
                        (
                            "Problem in dictionary molFractionFuel "
                            "(" + file_ + ")",
                            __FILE__,
                            __LINE__
                        );
                    }

                    if (modul)
                    {
                        data.insertCompositionFuelMol
                        (
                            tmp[element],
                            stod(tmp[element+1])
                        );
                    }
                }
            }
            else
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "Problem in dictionary molFractionFuel " 
                        "(" + file_ + ")",
                       __FILE__,
                        __LINE__
                    );
                }

                data.insertCompositionFuelMol
                (
                    tmp[0],
                    stod(tmp[1])
                );
            }
        }
    }
}


void AFC::PropertiesReader::massFractionFuel
(
    const stringList& fileContent,
    unsigned int line,
    Properties& data
)
{
    int dictBegin{-1};
    unsigned int dictEnd{0};
    
    findKeyword(dictBegin, dictEnd, fileContent, line);

    line = dictBegin+1;

    for(; line < dictEnd; line++)
    {
        //- Split string; delimiter ' ' 
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            //- If more elements in one line
            if (tmp.size() > 2)
            {
                forEach(tmp, element)
                {
                    unsigned int modul = element%2;

                    if (tmp[element+1].empty())
                    {
                        ErrorMsg
                        (
                            "Problem in dictionary massFractionFuel "
                            "(" + file_ + ")",
                            __FILE__,
                            __LINE__
                        );
                    }

                    if (modul)
                    {
                        data.insertCompositionFuelMass
                        (
                            tmp[element],
                            stod(tmp[element+1])
                        );
                    }
                }
            }
            else
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "Problem in dictionary massFractionFuel " 
                        "(" + file_ + ")",
                       __FILE__,
                        __LINE__
                    );
                }

                data.insertCompositionFuelMass
                (
                    tmp[0],
                    stod(tmp[1])
                );
            }
        }
    }
}


void AFC::PropertiesReader::control
(
    const stringList& fileContent,
    unsigned int line,
    Properties& data
)
{
    int dictBegin{-1};
    unsigned int dictEnd{0};
    
    findKeyword(dictBegin, dictEnd, fileContent, line);

    line = dictBegin+1;

    for(; line < dictEnd; line++)
    {
        //- Split string; delimiter ' ' 
        stringList tmp = splitStrAtWS(fileContent[line]);
        
        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            if (tmp[0] == "runTime")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No runTime specified in afcDict after keyword runTime "
                        "(" + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.insertRunTime(stod(tmp[1]));
            }
            else if (tmp[0] == "deltat")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No deltat specified in afcDict after keyword deltat "
                        "(" + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.insertDeltat(stod(tmp[1]));
            }
            else if (tmp[0] == "writeInterval")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No writeInterval specified in afcDict after keyword "
                        "writeInterval (" + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.insertWriteControlInterval(stod(tmp[1]));
            }
            else if (tmp[0] == "writeTime")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No writeTime specified in afcDict after keyword "
                        "writeTime (" + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.insertWriteControlTime(stod(tmp[1]));

            }
            else if (tmp[0] == "writeControl")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No writeControl specified in afcDict after keyword "
                        "writeControl (" + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                if
                (
                    tmp[1] != "endTime"
                 && tmp[1] != "interval"
                )
                {
                    Info<< "+ writeControl definition '" << tmp[1] << "' is not"
                        << " defined. Using endTime instead\n"
                        << endl;

                    data.insertWriteControl(tmp[1]);
                }
                else
                {
                    data.insertWriteControl(tmp[1]);
                }
            }
        }
    }
}


void AFC::PropertiesReader::input(const word inp, Properties& data)
{
    if (inp == "mol")
    {
        data.inputMol();
    }
    else if (inp == "mass")
    {
        data.inputMass();
    }
    else
    {
        ErrorMsg
        (
            "No fraction input found. Specify either mol or mass "
            "fraction in the afcDict.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::PropertiesReader::interpreter
(
    const stringList& fileContent,
    unsigned int line,
    Properties& data
)
{
    int dictBegin{-1};
    unsigned int dictEnd{0};
    
    findKeyword(dictBegin, dictEnd, fileContent, line);

    line = dictBegin+1;

    for(; line < dictEnd; line++)
    {
        //- Split string; delimiter ' ' 
        stringList tmp = splitStrAtWS(fileContent[line]);
        
        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            if (tmp[0] == "analyse")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No analyse keyword specified in afcDict ("
                        + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                //- Check input
                if (tmp[1] != "THERMO" && tmp[1] != "CHEMISTRY")
                {
                    ErrorMsg
                    (
                        "You can only use the keyword 'THERMO' "
                        "'CHEMSITRY' for analyse type",
                        __FILE__,
                        __LINE__
                    );
                }

                data.insertInterpreter(string(tmp[1]));
            }
        }
    }
}

// ************************************************************************* //
