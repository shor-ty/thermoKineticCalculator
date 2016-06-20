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
{
    private:


    public:

        //- Constructor
        TransportCalc();

        //- Destructor
        ~TransportCalc();

        // Reduced collision integrals

            //- Lenard 12-6 assumption, collision integral omega 2.2
            //  Neufeld et. al. (1972)
            scalar reducedCollisionIntegralOmega22
            (
                const scalar&
            );

            //- Lenard 12-6 assumption, collision integral omega 1.1
            //  Neufeld et. al. (1972)
            scalar reducedCollisionIntegralOmega11
            (
                const scalar&
            );

            //- Calculate density out of ideal gas law [kg/m^3]
            scalar rho
            (
                const scalar&,
                const scalar&,
                const scalar&
            );


        // Calculation functions for viscosity
        
            //- Kinetic calculation, pure species viscosity
            void viscosity
            (
                const Thermo&,
                TransportData&
            );

            //- Viscosity suggested by Neufeld et. al. (1972) [kg/m/s]
            scalar viscosityNeufeld
            (
                const word&,
                const scalar&,
                const Thermo&,
                const TransportData&
            );

            //- Viscosity suggested by Chung
            void viscosityChung
            (
                const word&
            );


        // Calculation functions for thermal conductivity

            //- Kinetic calculation, pure species thermal conductivity
            void thermalConductivity
            (
                const Thermo&,
                TransportData&
            );

            //- Thermal conductivity suggested by Warnatz [W/m/K]
            scalar thermalConductivityWarnatz3
            (
                const word&,
                const scalar&,
                const Thermo&,
                const TransportData&
            );

            //- Thermal conductivity suggested by Warnatz [W/m/K]
            scalar thermalConductivityWarnatz
            (
                const word&,
                const scalar&,
                const Thermo&,
                const TransportData&
            );


        // Calculation functions for binary diffusivity

            //- Kinetic calculation, pure species binary diffusivity
            void binaryDiffusivity
            (
                const Thermo&,
                TransportData&
            );



        // Return functions

            //- Return species as wordList
            wordList species() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // TransportCalc_hpp included

// ************************************************************************* //
