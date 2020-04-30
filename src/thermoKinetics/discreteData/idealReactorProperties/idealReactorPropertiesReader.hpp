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
    TKC::IdealReactorPropertiesReader

Description
    Reading the TKCDict file

SourceFiles
    propertiesReader.cpp

\*---------------------------------------------------------------------------*/

#ifndef IdealReactorPropertiesReader_hpp
#define IdealReactorPropertiesReader_hpp

#include "stringManipulator.hpp"
#include "idealReactorProperties.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                      Class IdealReactorPropertiesReader Declaration
\*---------------------------------------------------------------------------*/

class IdealReactorPropertiesReader
:
    public StringManipulator
{
    private:

        // Private data

            //- File name of ideal reactor parameter/setup
            string file_{"idealHomogeneousReactorDict"};

    public:

        // Constructor and Destructor

            //- Constructor with file path
            IdealReactorPropertiesReader(const string);

            //- Destructor
            ~IdealReactorPropertiesReader();


        // Member functions

            //- Read properties file and delegate data
            void read(IdealReactorProperties&);


        // Helper functions

            //- Find line number of keyword
            void findKeyword
            (
                int&,
                unsigned int&,
                const stringList&,
                unsigned int&
            );


        // Data manipulation functions

            //- Reading concentration dictionary
            void initialSpeciesData
            (
                const stringList&,
                unsigned int,
                IdealReactorProperties&
            );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
