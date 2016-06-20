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
    AFC::Properties
    
Description
    Abstract AFC::Properties class for properties data and calculation

SourceFiles
    properties.cpp

\*---------------------------------------------------------------------------*/

#ifndef Properties_hpp
#define Properties_hpp

#include "propertiesReader.hpp"
#include "thermo.hpp"
#include "chemistry.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Properties Declaration
\*---------------------------------------------------------------------------*/

class Properties
{
    private:


        // Debug
        const bool debug_{false};

        // Private Data
        
            //- Mixture fraction discrete points
            int nZPoints_{0};

            //- Varianz of mixture fraction discrete points
            int nZvarPoints_{0};

            //- Enthalpy defects [J/kg]
            scalarField defects_;

            //- ScalarDissipation rate [1/s]
            scalarField sDRs_;

            //- Oxidizer species
            wordList speciesOxidizer_;

            //- Composition of oxidizer mol fraction X [-]
            map<word, scalar> oxidizerX_;

            //- Composition of oxidizer mass fraction Y [-]
            map<word, scalar> oxidizerY_;

            //- Fuel species
            wordList speciesFuel_;

            //- Composition of fuel mol fraction X [-]
            map<word, scalar> fuelX_;

            //- Composition of fuel mass fraction Y [-]
            map<word, scalar> fuelY_;

            //- Temperature of oxidizer stream [K]
            scalar TOxidizer_{0};

            //- Temperature of fuel stream [K]
            scalar TFuel_{0};

            //- Pressure at which the calculation take place [Pa]
            scalar p_{0};


        // Boolean

            //- Input uses mol or mass fraction
            bool inputMol_{false};

            bool inputMass_{false};

            
        // Algorithm settings

            //- End time of calculation
            //  If time is high -> ddt(rho, Z)  == 0 steady state
            //  If time is low  -> ddt(rho, Z)  != 0 transient
            scalar runTime_{1e5};

            //- Write control of solution
            //  
            //      + endTime
            //      + iteration
            word writeControl_{"endTime"};

            //- Write control interval
            scalar writeControlInterval_{0};

            //- DeltaT
            scalar deltaT_{0};


        // Class object
            
            //- Reference to thermo object (XtoY and YtoX)
            const Thermo& thermo_;

            //- Reference to chemistry object (XtoY and YtoX)
            const Chemistry& chemistry_;


    public:

        //- Constructor
        Properties
        (
            const string&,
            Thermo&,
            const Chemistry&
        );

        //- Destructor
        ~Properties();


        // Member functions
        
            //- Insert nZPoints_
            void insertMFPoints
            (
                const int&
            );
        
            //- Insert nZvarPoints_
            void insertVMFPoints
            (
                const int&
            );

            //- Insert 
            void insertEnthalpyDefects
            (
                const scalar&
            );

            //- Insert scalar dissipation rates
            void insertScalarDissipationRates
            (
                const scalar&
            );

            //- Insert temperature oxidizer
            void insertTemperatureOxidizer
            (
                const scalar&
            );

            //- Insert temperature fuel
            void insertTemperatureFuel
            (
                const scalar&
            );
            
            //- Insert oxidizer mol composition 
            void insertCompositionOxidizerMol
            (
                const word&,
                const scalar&
            );
            
            //- Insert oxidizer mass composition 
            void insertCompositionOxidizerMass
            (
                const word&,
                const scalar&
            );
            
            //- Insert fuel mol composition 
            void insertCompositionFuelMol
            (
                const word&,
                const scalar&
            );
            
            //- Insert fuel mass composition 
            void insertCompositionFuelMass
            (
                const word&,
                const scalar&
            );

            //- Insert endTime of calculation
            void insertRunTime
            (
                const scalar&
            );

            //- Insert control of write the solution
            void insertWriteControl
            (
                const word&
            );
            
            //- Insert write control (save each x iterations)
            void insertWriteControlInterval
            (
                const scalar&
            );

            //- Insert deltaT
            void insertDeltaT
            (
                const scalar&
            );

            //- Insert pressure
            void insertPressure
            (
                const scalar&
            );

            //- Insert bool for inputMol
            void inputMol();

            //- Insert bool for inputMass
            void inputMass();


        // Other functions

            //- Check if all data are set
            void check();

            //- Mol fraction to mass fraction
            void XtoY();

            //- Mass fraction to mol fraction
            //  + word -> oxidizer (O) or fuel (F)
            /*void YtoX
            (
                const word 
            );*/


        // Return functions

            //- Return oxidizer species
            wordList speciesOxidizer() const;

            //- Return oxidizer mol fraction
            map<word, scalar> oxidizerCompMol() const;
            //
            //- Return oxidizer mass fraction
            map<word, scalar> oxidizerCompMass() const;

            //- Return fuel species
            wordList speciesFuel() const;

            //- Return fuel mol fraction
            map<word, scalar> fuelCompMol() const;

            //- Return fuel mass fraction
            map<word, scalar> fuelCompMass() const;

            //- Return scalar dissipation rates [1/s]
            scalarField sDRs() const;

            //- Return enthalpy defects [J/kg]
            scalarField defects() const;

            //- Return the number of discrete points of Z
            int nZPoints() const;

            //- Return the number of discrete points of Zvar
            int nZvarPoints() const;

            //- Return oxidizer temperature [K]
            scalar oxidizerTemperature() const;

            //- Return fuel temperature [K]
            scalar fuelTemperature() const;

            //- Return runTime [s]
            scalar runTime() const;

            //- Return deltaT
            scalar deltaT() const;

            //- Return number of defects
            unsigned int nDefects() const;

            //- Return defect value [J/kg]
            scalar defect
            (
                const int&
            ) const;

            //- Return pressure [Pa]
            scalar p() const;

            //- Return word of input (mol or mass)
            word input() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Properties_hpp included

// ************************************************************************* //
