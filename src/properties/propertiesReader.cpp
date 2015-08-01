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

AFC::PropertiesReader::PropertiesReader
(
    const string& file
)
:
    file_(file)

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::PropertiesReader::~PropertiesReader()
{
//    if (debug)
//    {
        Info<< "Destructor PropertiesReader\n" << endl;
//    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AFC::PropertiesReader::read
(
    Properties& data
)
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
            if (tmp[0] == "mixtureFractionPoints")
            {
                if (tmp[1].empty())
                {
                    FatalError
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
                    FatalError
                    (
                        "No varianzOfMixtureFractionPoints specified ("
                        + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.insertVMFPoints(stoi(tmp[1]));
            }
            else if (tmp[0] == "scalarDissipationRates")
            {
                scalarDissipationRates(fileContent, line, data);
            }
            else if (tmp[0] == "temperatureOxidizer")
            {
                data.insertTemperatureOxidizer(stod(tmp[1]));
            }
            else if (tmp[0] == "temperatureFuel")
            {
                data.insertTemperatureFuel(stod(tmp[1]));
            }
            else if (tmp[0] == "molFractionOxidizer")
            {
                molFractionOxidizer(fileContent, line, data);
            }
            else if (tmp[0] == "molFractionFuel")
            {
                molFractionFuel(fileContent, line, data);
            }
        }
    }

    data.check();
}


// * * * * * * * * * * * * * * Helper functions  * * * * * * * * * * * * * * //

void AFC::PropertiesReader::findKeyword
(
    int& start,
    unsigned int& end,
    const stringList& fileContent,
    unsigned int line
)
{
    for (; line < fileContent.size(); line++)
    {
        //- Split string; delimiter ' ' 
        stringList tmp = splitStrAtWS(fileContent[line]);

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


// * * * * * * * * * * * * Data manipulation functions * * * * * * * * * * * //

void AFC::PropertiesReader::scalarDissipationRates
(
    const stringList& fileContent,
    unsigned int& line,
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
                forAll(tmp, element)
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
    unsigned int& line,
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
                forAll(tmp, element)
                {
                    unsigned int modul = element%2;

                    if (modul)
                    {
                        if (tmp[element+1].empty())
                        {
                            FatalError
                            (
                                "    Error in dictionary molFractionOxidizer "
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
                    FatalError
                    (
                        "    Error in dictionary molFractionOxidizer "
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


void AFC::PropertiesReader::molFractionFuel
(
    const stringList& fileContent,
    unsigned int& line,
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
                forAll(tmp, element)
                {
                    unsigned int modul = element%2;

                    if (tmp[element+1].empty())
                    {
                        FatalError
                        (
                            "    Error in dictionary molFractionFuel "
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
                    FatalError
                    (
                        "    Error in dictionary molFractionFuel " 
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


// ************************************************************************* //
