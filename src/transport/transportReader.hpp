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
    AFC::TransportReader    

Description
    Reading the chemkin III file

SourceFiles
    transportReader.cpp

\*---------------------------------------------------------------------------*/

#ifndef TransportReader_hpp
#define TransportReader_hpp

#include "stringManipulator.hpp"
#include "transportData.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                      Class TransportReader Declaration
\*---------------------------------------------------------------------------*/

class TransportReader
:
    public StringManipulator
{
    private:

        // Private data

            //- Transport file
            string file_;


    public:

        // Constructor and Destructor

            //- Constructor with file string and Transport:: obj adress
            TransportReader
            (
                const string&
            );

            //- Destructor
            ~TransportReader();

        
        // Member functions

            //- Read transport file and return pointer to TransportData:: obj
            void read
            (
                TransportData& 
            );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
