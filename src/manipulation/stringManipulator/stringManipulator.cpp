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

#include "stringManipulator.hpp"
#include <fstream>
#include <iterator>
#include <sstream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::StringManipulator::StringManipulator()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::StringManipulator::~StringManipulator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const TKC::stringList TKC::StringManipulator::readFile(const string filePath)
{
    std::ifstream file;

    file.open(filePath.c_str(), std::ios::in);

    if (!file.good())
    {
        ErrorMsg
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


const TKC::stringList TKC::StringManipulator::splitStrAtWS(const string str)
{
    //- Split string at whitespace
    std::istringstream tmp(str);

    stringList sF
    {
        std::istream_iterator<std::string>{tmp},
        std::istream_iterator<std::string>{}
    };

    return sF;
}


const TKC::stringList TKC::StringManipulator::splitStrAtDelimiter
(
    const string str,
    const char delimiter
)
{
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


const TKC::string
TKC::StringManipulator::removeAtEnd(const string str, const string str2)
{
    //- Find position of str2
    size_t pos = str.find(str2);

    //- Return string von pos 0 to pos
    return std::move(str.substr(0, pos));
}


void TKC::StringManipulator::removeFirstChar(string& str)
{
    str.replace(0, 1, "");
}


void TKC::StringManipulator::removeLastChar(string& str)
{
    str.replace(str.size()-1, 1, "");
}


void TKC::StringManipulator::removeComment(string& str, const string comment)
{
    //- Search comment position
    if (!str.empty())
    {
        size_t pos = str.find(comment);

        str = str.substr(0, pos);
    }
}


const TKC::string TKC::StringManipulator::extractBetweenKeys
(
    const string& str,
    const string key1,
    const string key2
)
{
    std::size_t pos1 = str.find(key1);
    std::size_t pos2 = str.find(key2);

    return std::move(str.substr(pos1, pos2-pos1+1));
}


// ************************************************************************* //
