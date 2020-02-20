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

        // Warnings
        bool warnings{false};


    public:

        //- Constructor
        ThermoCalc();

        //- Destructor
        ~ThermoCalc();


        // Calculation functions

            //- Calculate mean molecular weight out of the mol fraction
            //  [kg/mol]
            scalar MWmeanX
            (
                const map<word, scalar>&,
                const map<word, scalar>&
            ) const;
            
            //- Calculate mean molecular weight out of the mass fraction
            //  [kg/mol]
            scalar MWmeanY
            (
                const map<word, scalar>&,
                const map<word, scalar>&
            ) const;
            
            //- Calculate mean molecular weight out of the concentration
            //  [kg/mol]
            scalar MWmeanC
            (
                const map<word, scalar>&,
                const map<word, scalar>&
            ) const;

            //- Calculate rho for species based on the ideal gas law
            //  [kg/m^3]
            scalar rho(const scalar, const scalar, const scalar) const;

            //- Calculate mean density based on the mean molecular weight
            //  [kg/m^3]
            scalar rhoMean(const scalar, const scalar, const scalar) const;

            //- Calculate complete concentration C [mol/m^3]
            scalar C(const scalar, const scalar) const;
            

            //- Calculate temperature depended specific heat capacity via
            //  NASA polynomials of species s for constant pressure
            //  [J/mol/K] 
            scalar cp(const word, const scalar, const ThermoData&) const;

            //- Calculate temperature depended specific heat capcity via
            //  NASA polynomials of species s for constant volume
            //  [J/mol/K] 
            scalar cv(const word, const scalar, const ThermoData&) const;

            //- Calculate temperature dependend enthalpy via 
            //  NASA polynomials of species s
            //  [J/mol]
            scalar H(const word, const scalar, const ThermoData&) const;

            //- Calculate the enthalpy difference from H and Hf via
            //  NASA polynomials for species s
            //  [J/mol]
            scalar dHf(const word, const scalar, const ThermoData&) const;

            //- Calculate temperature dependend entropy via
            //  NASA polynomials for species s
            //  [J/mol/K] 
            scalar S(const word, const scalar, const ThermoData&) const;

            //- Calculate temperature dependend free GIBBS energy via
            //  NASA polynomials for species s
            //  [J/mol/K] 
            scalar G(const word, const scalar, const ThermoData&) const;

            //- Calculate the free GIBBS energy difference from G and G(298)
            //  [J/mol]
            scalar dGf(const word, const scalar, const ThermoData&) const;

            //- Calculate temperature dependend free GIBBS energy via
            //  NASA polynomials of the mixture (mean value)
            //  [J/mol] 
            scalar G(const scalar, const scalar, const scalar) const;

            //- Calculate formation enthalpy Hf == H(298) of species s
            //  [J/mol]
            scalar Hf(const word, const ThermoData&) const;

            //- Calculate formation free GIBBS energy Gf == G(298) of species s
            //  [J/mol]
            scalar Gf(const word, const ThermoData&) const;


        // Return functions
        
            //- Get correct NASA coeffs
            scalarField getCoeffs
            (
                const word,
                const scalar,
                const ThermoData&
            ) const;

            //- Checking for correct temperature range
            bool whichTempRange
            (
                const word,
                const scalar,
                const ThermoData&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ThermoCalc_hpp included

// ************************************************************************* //
