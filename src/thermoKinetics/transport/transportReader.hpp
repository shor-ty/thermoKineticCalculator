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
    TKC::TransportReader    

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

namespace TKC
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

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
