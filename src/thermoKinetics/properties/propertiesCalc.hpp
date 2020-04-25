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
    TKC::PropertiesCalc
    
Description
    Abstract TKC::PropertiesCalc class for manipulation, calculation 

SourceFiles
    propertiesCalc.cpp

\*---------------------------------------------------------------------------*/

#ifndef PropertiesCalc_hpp
#define PropertiesCalc_hpp

#include "propertiesData.hpp"
#include "thermo.hpp"
#include "chemistry.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                            Class Properties Declaration
\*---------------------------------------------------------------------------*/

class PropertiesCalc
:
    public PropertiesData
{
    private:

        // Referecne to class object
            
            //- Reference to thermo object (XtoY and YtoX)
            const Thermo& thermo_;

            //- Reference to chemistry object (XtoY and YtoX)
            const Chemistry& chemistry_;
        

    public: 

        //- Constructor
        PropertiesCalc(const string, const Thermo&, const Chemistry&);

        //- Destructor
        ~PropertiesCalc();



        // Calculation functions
        
        /*
            //- Calculate the minimum oxigen needed at stochiometric condition
            static scalar Omin
            (
                const map<word, scalar>&,
                const Thermo&
            );
        
            //- Calculate the stochiometric mixture fraction Zst
            static scalar Zst
            (
                const map<word, scalar>&,
                const scalar,
                const scalar,
                const Thermo&
            );

            //- Calculate the adiabatic flame temperature and return it [K]
            static scalar adiabateFlameTemperature
            (
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
            static map<word, scalar> elementMassFraction
            (
                const map<word, map<word, unsigned int>>&,
                const map<word, scalar>&,
                const Thermo&
            );

            //- Create the element mole fraction of the given chemical
            //  species; used for oxidizer and fuel site
            static map<word, scalar> elementMolFraction
            (
                const map<word, scalar>&,
                const Thermo&
            );
            */
              
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // PropertiesCalc_hpp included

// ************************************************************************* //
