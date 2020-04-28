/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     T hermo       | Open Source Thermo-Kinetic Library
  c-o-o-c     K iknetic     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Constructor.

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

\*---------------------------------------------------------------------------*/

#include "idealReactorPropertiesReader.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::IdealReactorPropertiesReader::IdealReactorPropertiesReader
(
    const string file
)
{
    if (!file.empty())
    {
        file_ = file;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::IdealReactorPropertiesReader::~IdealReactorPropertiesReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void TKC::IdealReactorPropertiesReader::read(IdealReactorProperties& data)
{
    Info<< " c-o Reading dictonary\n"
        << "     >> " << file_ << "\n" << endl;

    const auto fileContent = readFile(file_);

    //- Reading afcDict
    for (unsigned int line=0; line < fileContent.size(); line++)
    {
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if (!tmp.empty() && tmp[0][0] != '!')
        {
            if (tmp[0] == "inertGas")
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
            else if (tmp[0] == "pressure")
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

                data.p(stod(tmp[1]));
            }
            else if (tmp[0] == "temperature")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No temperature specified (" + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                data.T(stod(tmp[1]));
            }
            else if (tmp[0] == "concentration")
            {
                data.inputMode("concentration");
                initialSpeciesData(fileContent, line, data);
            }
            else if (tmp[0] == "moleFraction")
            {
                data.inputMode("mole");
                initialSpeciesData(fileContent, line, data);
            }
            else if (tmp[0] == "massFraction")
            {
                data.inputMode("mass");
                initialSpeciesData(fileContent, line, data);
            }
            else if (tmp[0] == "thermodynamic")
            {
                data.thermo(tmp[1]);
            }
            else if (tmp[0] == "chemistry")
            {
                data.chemistry(tmp[1]);
            }
            else if (tmp[0] == "transport")
            {
                data.transport(tmp[1]);
            }
            else if (tmp[0] == "interprete")
            {
                if (tmp[1] == "true" || tmp[1] == "yes")
                {
                    data.interprete(true);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * Helper functions  * * * * * * * * * * * * * * //

void TKC::IdealReactorPropertiesReader::findKeyword
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
        if (!tmp.empty() && tmp[0][0] != '!')
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
            if (start != 0 && end != 0)
            {
                break;
            }
        }
    }
}


// * * * * * * * * * * * * Data manipulation functions * * * * * * * * * * * //


void TKC::IdealReactorPropertiesReader::initialSpeciesData
(
    const stringList& fileContent,
    unsigned int line,
    IdealReactorProperties& data
)
{
    int dictBegin{-1};
    unsigned int dictEnd{0};

    findKeyword(dictBegin, dictEnd, fileContent, line);

    line = dictBegin+1;

    for(; line < dictEnd; line++)
    {
        //- Line content
        string lineContent = fileContent[line];

        //- Remove any comments '!'
        removeComment(lineContent);

        //- Split string; delimiter ' '
        stringList tmp = splitStrAtWS(lineContent);

        //- If line is not empty and no comment, proceed
        if (!tmp.empty())
        {
            //- If more elements in one line
            if (tmp.size() > 2)
            {
                //- The modulus has to be equal to zero
                if (!tmp.size() % 2)
                {
                    ErrorMsg
                    (
                        "Problem in dictionary ",
                        __FILE__,
                        __LINE__
                    );
                }

                forEach(tmp, element)
                {
                    unsigned int modul = element%2;

                    Info<< "element: " << element << endl;

                    /*
                    if (modul)
                    {
                        data.insertCompositionOxidizerMass
                        (
                            tmp[element],
                            stod(tmp[element+1])
                        );
                    }
                    */
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
            }
        }
    }
}


// ************************************************************************* //
