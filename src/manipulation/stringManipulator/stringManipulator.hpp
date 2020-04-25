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

Class
    TKC::StringManipulator

Description
    Class for string manipulations and file operations

SourceFiles
    stringManipulation.cpp

\*---------------------------------------------------------------------------*/

#ifndef StringManipulator_hpp
#define StringManipulator_hpp

#include "definitions.hpp"
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                      Class StringManipulator Declaration
\*---------------------------------------------------------------------------*/

class StringManipulator
{
    public:

        //- Constructor
        StringManipulator();

        //- Destructor
        ~StringManipulator();


        // Member function

            //- Open file and return file content as stringList
            const stringList readFile(const string);

            //- Split string; delimiter is space ' '
            const stringList splitStrAtWS (const string);

            //- Split string; delimiter is a character
            const stringList splitStrAtDelimiter(const string, const char);

            //- Split string; delimiter is a string
            const stringList splitStrAtDelimiter(const string, const string);

            //- Remove str2 from str and return str from 0 - pos
            const string removeAtEnd(const string, const string);

            //- Remove first character from string
            void removeFirstChar(string&);

            //- Remove last character from string
            void removeLastChar(string&);

            //- Remove the comment in the string
            void removeComment(string&, const string comment = "!");

            //- Extract string between two keys
            const string extractBetweenKeys
            (
                const string&,
                const string key1 = "(",
                const string key2 = ")"
            );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // StringManipulator_hpp included

// ************************************************************************* //
