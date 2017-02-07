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

#include "stringManipulator.hpp"
#include <fstream>
#include <iterator>
#include <sstream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::StringManipulator::StringManipulator()
{
    if (debug_)
    {
        Info<< "Constructor StringManipulator\n" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::StringManipulator::~StringManipulator()
{
    if (debug_)
    {
        Info<< "Destructur StringManipulator\n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const AFC::stringList AFC::StringManipulator::readFile
(
    const string& filePath
)
{
    if (debug_)
    {
        Info<< " --> AFC::StringManipulator::readFile" << endl;
    }

    std::ifstream file;

    file.open(filePath.c_str(), std::ios::in);

    if (!file.good())
    {
        FatalError
        (
            "    File \"" + filePath + "\" can not be opened.",
            __FILE__,
            __LINE__
        );
    }

    string lineContent;
    stringList fileContent;

    //- Store all lines in fileContent and return
    while(!file.eof())
    {
        std::getline(file, lineContent);
        fileContent.push_back(lineContent);
    }

    //- Remove last entry
    fileContent.resize(fileContent.size()-1);

    //- Close file
    file.close();

    return fileContent;
}


const AFC::stringList AFC::StringManipulator::splitStrAtWS
(
    const string& str
)
{
    if (debug_)
    {
//        Info<< " --> AFC::StringManipulator::splitStrAtWS" << endl;
    }

    //- Split string at whitespace
    std::istringstream tmp(str);

    stringList sF
    {
        std::istream_iterator<std::string>{tmp},
        std::istream_iterator<std::string>{}
    };

    return sF;
}


const AFC::stringList AFC::StringManipulator::splitStrAtDelimiter
(
    const string& str,
    const char& delimiter
)
{
    if (debug_)
    {
        Info<< " --> AFC::StringManipulator::splitStrAtDelimiter" << endl;
        Info<< "String to split: " << str << endl;
        Info<< "Delimiter is: " << delimiter << endl;
    }

    //- Split string at delimiter and return the stringList
    std::stringstream tmp(str);
    string element;
    stringList elements;

    while (std::getline(tmp, element, delimiter))
    {
        elements.push_back(element);
    }

    return elements;
}


const AFC::string AFC::StringManipulator::removeAtEnd
(
    const string& str,
    const string& str2
)
{
    if (debug_)
    {
        Info<< " --> AFC::StringManipulator::removeAtEnd" << endl;
    }

    //- Find position of str2
    std::size_t pos = str.find(str2);

    //- Return string von pos 0 to pos
    return str.substr(0,pos);
}


void AFC::StringManipulator::removeFirstChar
(
    string& str
)
{
    str.replace(0,1,"");
}


// ************************************************************************* //
