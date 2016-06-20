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

Class
    AFC::StringManipulator

Description
    Class for string manipulations and file operations 

SourceFiles
    stringManipulation.cpp

\*---------------------------------------------------------------------------*/

#ifndef StringManipulator_hpp
#define StringManipulator_hpp

#include "typedef.hpp"
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                      Class StringManipulator Declaration
\*---------------------------------------------------------------------------*/

class StringManipulator
{
    protected:

        // Debug switch
        bool debug_{false};


    private:


    public:

        //- Constructor
        StringManipulator();

        //- Destructor
        ~StringManipulator();


        // Member function

            //- Open file and return file content as stringList
            const stringList readFile
            (
                const string& 
            );

            //- Split string; delimiter ' '
            const stringList splitStrAtWS 
            (
                const string&
            );

            //- Split string; delimiter as parameter
            const stringList splitStrAtDelimiter
            (
                const string&,
                const char&
            );

            //- Split string; delimiter as parameter
            const stringList splitStrAtDelimiter
            (
                const string&,
                const string&
            );

            //- Remove str2 from str and return str from 0 - pos
            const string removeAtEnd
            (
                const string&,
                const string&
            );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // StringManipulator_hpp included

// ************************************************************************* //
