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
    TKC::PropertiesReader

Description
    Reading the TKCDict file

SourceFiles
    propertiesReader.cpp

\*---------------------------------------------------------------------------*/

#ifndef PropertiesReader_hpp
#define PropertiesReader_hpp

#include "stringManipulator.hpp"
#include "propertiesData.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                      Class PropertiesReader Declaration
\*---------------------------------------------------------------------------*/

class PropertiesReader
:
    public StringManipulator
{
    private:

        // Private data

            //- Properties file
            const string file_;


    public:

        // Constructor and Destructor

            //- Constructor with file string and Properties:: obj adress
            PropertiesReader(const string);

            //- Destructor
            ~PropertiesReader();


        // Member functions

            //- Read properties file and delegate data
            void read(PropertiesData&);

            //- Return the path to the thermo, transport and chemistry file 
            string path(const string);
            
            //- Return if the data should be interpreted
            bool interprete();

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

            //- Reading enthalpy defects dictionary
            void enthalpyDefects
            (
                const stringList&,
                unsigned int,
                PropertiesData&
            );

            //- Reading scalar dissipation rate dictionary
            void scalarDissipationRates
            (
                const stringList&,
                unsigned int,
                PropertiesData&
            );

            //- Reading mol fraction oxidizer dictionary
            void molFractionOxidizer
            (
                const stringList&,
                unsigned int,
                PropertiesData&
            );

            //- Reading mass fraction oxidizer dictionary
            void massFractionOxidizer
            (
                const stringList&,
                unsigned int,
                PropertiesData&
            );

            //- Reading mol fraction fuel dictionary
            void molFractionFuel
            (
                const stringList&,
                unsigned int,
                PropertiesData&
            );

            //- Reading mass fraction fuel dictionary
            void massFractionFuel
            (
                const stringList&,
                unsigned int,
                PropertiesData&
            );

            //- Reading control dictionary
            void control(const stringList&, unsigned int, PropertiesData&);

            //- Set wheater input is mol or mass fraction
            void input(const word, PropertiesData&);

            //- Read interpreter dictionary
            void interpreter(const stringList&, unsigned int, PropertiesData&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
