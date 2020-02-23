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
    AFC::PropertiesReader

Description
    Reading the AFCDict file

SourceFiles
    propertiesReader.cpp

\*---------------------------------------------------------------------------*/

#ifndef PropertiesReader_hpp
#define PropertiesReader_hpp

#include "stringManipulator.hpp"
#include "properties.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{


class Properties;

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
            string file_;


    public:

        // Constructor and Destructor

            //- Constructor with file string and Properties:: obj adress
            PropertiesReader(const string);

            //- Destructor
            ~PropertiesReader();


        // Member functions

            //- Read properties file and delegate data
            void read(Properties&);

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
            void enthalpyDefects(const stringList&, unsigned int, Properties&);

            //- Reading scalar dissipation rate dictionary
            void scalarDissipationRates
            (
                const stringList&,
                unsigned int,
                Properties&
            );

            //- Reading mol fraction oxidizer dictionary
            void molFractionOxidizer
            (
                const stringList&,
                unsigned int,
                Properties&
            );

            //- Reading mass fraction oxidizer dictionary
            void massFractionOxidizer
            (
                const stringList&,
                unsigned int,
                Properties&
            );

            //- Reading mol fraction fuel dictionary
            void molFractionFuel(const stringList&, unsigned int, Properties&);

            //- Reading mass fraction fuel dictionary
            void massFractionFuel(const stringList&, unsigned int, Properties&);

            //- Reading control dictionary
            void control(const stringList&, unsigned int, Properties&);

            //- Set wheater input is mol or mass fraction
            void input(const word, Properties&);

            //- Read interpreter dictionary
            void interpreter(const stringList&, unsigned int, Properties&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
