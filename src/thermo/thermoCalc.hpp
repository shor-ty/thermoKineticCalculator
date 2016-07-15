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
    AFC::ThermoCalc class for thermo calculation. This class provides all 
    thermo calculations that are necessary in the AFC project. All calculation
    functions are given below (ideal gas assumption):

    \f[ \bar{M} = X_i M_i = \left(\frac{Y_i}{M_i}\right)^{-1} \f]

    \f[ \rho = \frac{p \bar{M}}{RT} \f]

    \f[ \rho = \sum_s^n\left(\frac{X_s p}{R M_s \sum_r^n X_r} \right)\f]

    \f[ \rho = \frac{p}{RT \sum_{i=1}^n \frac{Y_i}{M_i}} \f]
     
    \f[ c = \frac{p}{RT} \f]

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

            //- Calculate mean molecular weight out of the mol fraction [g/mol]
            scalar MmeanX
            (
                const map<word, scalar>&,
                const map<word, scalar>&
            ) const;
            
            //- Calculate mean molecular weight out of the mass fraction [g/mol]
            scalar MmeanY
            (
                const map<word, scalar>&,
                const map<word, scalar>&
            ) const;
            
            //- Calculate mean molecular weight out of the concentration [g/mol]
            scalar MmeanC
            (
                const map<word, scalar>&,
                const map<word, scalar>&
            ) const;

            //- Calculate mean density based on the mean molecular weight [g/m^3]
            scalar rhoMean
            (
                const scalar&,
                const scalar&,
                const scalar&
            ) const;

            //- Calculate complete concentration C [mol/m^3]
            scalar C
            (
                const scalar&,
                const scalar&
            ) const;
            

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

            //- Calculate the enthalpy difference from H and Hf (NASA) [J/mol]
            scalar dH
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

            //- Calculate the free GIBBS energy difference from G and G(298)
            //  [J/mol]
            scalar dG
            (
                const word&,
                const scalar&,
                const ThermoData&
            ) const;


            //- Calculate temperature dependend free GIBBS energy (NASA)
            //  of mixture (mean) [J/mol] 
            scalar G
            (
                const scalar&,
                const scalar&,
                const scalar&
            ) const;

            //- Calculate formation enthalpy Hf == H(298) [J/mol]
            scalar Hf
            (
                const word&,
                const ThermoData&
            ) const;

            //- Calculate formation free GIBBS energy Gf == G(298) [J/mol]
            scalar Gf
            (
                const word&,
                const ThermoData&
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
