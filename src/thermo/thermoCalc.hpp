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
    AFC::ThermoCalc
    
Description
    Abstract AFC::ThermoCalc class for thermo calculation

SourceFiles
    thermo.cpp

\*---------------------------------------------------------------------------*/

#ifndef ThermoCalc_hpp
#define ThermoCalc_hpp

#include "thermoData.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class ThermoCalc Declaration
\*---------------------------------------------------------------------------*/

class ThermoCalc
{
    private:

        // Debug
        bool debug{true};

        // Warnings
        bool warnings{false};


    public:

        //- Constructor
        ThermoCalc();

        //- Destructor
        ~ThermoCalc();


        // Calculation functions

            //- Calculate temperature dependend heat capacity (NASA) [J/mol/K]
            //  of species s
            scalar cp
            (
                const word&,
                const scalar&,
                const ThermoData&
            ) const;

            //- Calculate temperature dependend enthalpy (NASA) [J/mol]
            //  of species s
            scalar H
            (
                const word&,
                const scalar&,
                const ThermoData&
            ) const;

            //- Calculate temperature dependend entropy (NASA) [J/mol/K] 
            //  of species s
            scalar S
            (
                const word&,
                const scalar&,
                const ThermoData&
            ) const;

            //- Calculate temperature dependend free GIBBS energy (NASA)
            //  of species s [J/mol/K] 
            scalar G
            (
                const word&,
                const scalar&,
                const ThermoData&
            ) const;

            //- Calculate temperature dependend free GIBBS energy (NASA)
            //  of mixture (mean) [J/mol/K] 
            scalar G
            (
                const scalar&,
                const scalar&,
                const scalar&
            ) const;


        // Return functions
        
            //- Get correct NASA coeffs
            scalarField getCoeffs
            (
                const word&,
                const scalar&,
                const ThermoData&
            ) const;

            //- Return true or false to choose the correct NASA polys
            //  of species s
            bool whichTempRange
            (
                const word&,
                const scalar&,
                const ThermoData&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ThermoCalc_hpp included

// ************************************************************************* //
