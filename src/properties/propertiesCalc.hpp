/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
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
    AFC::Properties
    
Description
    Abstract AFC::Properties class for propertiesCalc data and calculation

SourceFiles
    propertiesCalc.cpp

\*---------------------------------------------------------------------------*/

#ifndef PropertiesCalc_hpp
#define PropertiesCalc_hpp

#include "properties.hpp"
#include "thermo.hpp"
#include "chemistry.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Properties Declaration
\*---------------------------------------------------------------------------*/

class PropertiesCalc
{
    private:

        // General private data
        

    public: 

        //- Constructor
        PropertiesCalc();

        //- Destructor
        ~PropertiesCalc();



        // Calculation functions
        
            //- Calculate the stochiometric mixture fraction Zst 
            static scalar Zst
            (
                const map<word, scalar>&,
                const map<word, scalar>&
            );

            //- Calculate the adiabatic flame temperature and return it [K]
            static scalar adiabateFlameTemperature
            (
                const Properties&,
                const Thermo&,
                const Chemistry&
            );

            //- Decompose the species and count the atoms
            static map<word, map<word, unsigned int>> elementDecomposition
            (
                const wordList&,
                const Thermo&
            );

            //- Create the element mass fraction of the given chemical
            //  species; used for oxidizer and fuel site
            static map<word, unsigned int> elementMassFraction
            (
                const map<word, unsigned int>&,
                const Thermo&
            );
              
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // PropertiesCalc_hpp included

// ************************************************************************* //
