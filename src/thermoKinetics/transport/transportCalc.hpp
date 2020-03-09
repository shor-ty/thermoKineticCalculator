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
    AFC::TransportCalcCalc
    
Description
    Abstract AFC::TransportCalcCalc class for transport calculation
    Mostly gas-kinetics are implemented and a fitting procedure    

SourceFiles
    transportCalc.cpp

\*---------------------------------------------------------------------------*/

#ifndef TransportCalc_hpp
#define TransportCalc_hpp

#include "transportData.hpp"
#include "thermo.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class TransportCalc Declaration
\*---------------------------------------------------------------------------*/

class TransportCalc
:
    public TransportData
{
    private:

        // Reference to thermo class
        const Thermo& thermo_;


    public:

        //- Constructor
        TransportCalc(const string, const Thermo&);

        //- Destructor
        ~TransportCalc();

        // Reduced collision integrals

            //- Lenard 12-6 assumption, collision integral omega 2.2
            //  Neufeld et. al. (1972)
            scalar reducedCollisionIntegralOmega22(const scalar) const;

            //- Lenard 12-6 assumption, collision integral omega 1.1
            //  Neufeld et. al. (1972)
            scalar reducedCollisionIntegralOmega11(const scalar) const;

            //- Calculate density out of ideal gas law [kg/m^3]
            scalar rho(const scalar, const scalar, const scalar) const;


        // Calculation functions for viscosity
        
            //- Kinetic calculation, pure species viscosity
            scalar viscosity
            (
                const word,
                const scalar,
                const word method = "Hirschfelder"
            ) const;

            //- Return the viscosity suggested by Hirschfelder et. al. (1954)
            //  [kg/m/s]
            scalar viscosityHirschfelder
            (
                const word,
                const scalar,
                const scalar,
                const scalar,
                const scalar
            ) const;

            //- Return the viscosity suggested by Chung
            //  [kg/m/s]
            scalar viscosityChung(const word, const scalar) const;

            //- Return the viscosity using the polynomial
            //  [kg/m/s]
            scalar viscosityPolynomial(const scalar, const scalarField&) const;

            //- Fitting function for viscosity
            void fitViscosity();


        // Calculation functions for thermal conductivity

            //- Kinetic calculation, pure species thermal conductivity
            scalar thermalConductivity
            (
                const word,
                const scalar,
                const word method = "Warnatz" 
            ) const;

            //- Return the thermal conductivity suggested by Warnatz in the
            //  book Verbrennung by Warnatz et. al.
            //  [W/m/K]
            scalar thermalConductivityWarnatz(const word, const scalar) const;

            //- Return the thermal conductivity suggested by Warnatz
            //  This formula is mentioned in the TRANSPORT summary of
            //  Chemkin Collection Release 3.6 (2000)
            //  [W/m/K]
            scalar thermalConductivityWarnatzCC(const word, const scalar) const;

            //- Return the thermal conductivity using the polynomial 
            scalar thermalConductivityPolynomial
            (
                const scalar,
                const scalarField&
            ) const;

            //- Fitting function for thermal conductivity
            void fitThermalConductivity();


        // Calculation functions for binary diffusivity

            //- Kinetic calculation, pure species binary diffusivity
            scalar binaryDiffusivity
            (
                const word,
                const word,
                const scalar,
                const word method = "ChapmanAndEnskog"
            ) const;

            //- Return the binary diffusivity suggested by Warnatz et. al.
            //  in the book Verbrennung by Warnatz et. al.
            scalar binaryDiffusivityWarnatz
            (
                const word,
                const word,
                const scalar
            ) const;

            //- Return the binary diffusivity suggested by Chapman and Enskog
            //  [m^2/s]
            scalar binaryDiffusivityChapmanAndEnskog
            (
                const word,
                const word,
                const scalar
            ) const;

            //- Return the binary diffusivity suggested in the TRANSPORT
            //  summary of the Chemkin Collection Release 3.6 (2000)
            scalar binaryDiffusivityCC
            (
                const word,
                const word,
                const scalar
            ) const;

            //- Return the binary diffusivity using the polynomial
            //  [m^2/s]
            scalar binaryDiffusivityPolynomial
            (
                const scalar,
                const scalarField&
            ) const;

            //- Fitting function for binary diffusivity 
            void fitBinaryDiffusivity();


        // Additional Functions

            //- Return the temperature dependence of ZRot
            scalar F(const scalar, const scalar) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // TransportCalc_hpp included

// ************************************************************************* //
