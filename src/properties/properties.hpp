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

        // Private Data
        
            //- Mixture fraction discrete points
            int mfPoints_{0};

            //- Varianz of mixture fraction discrete points
            int vmfPoints_{0};

            //- Enthalpy defects
            scalarField defects_;

            //- ScalarDissipation rate
            scalarField sDRs_;

            //- Oxidizer species
            wordList speciesOxidizer_;

            //- Composition of oxidizer [mol] mol fraction X
            map<word, scalar> oxidizerX_;

            //- Composition of oxidizer [kg] mass fraction Y
            map<word, scalar> oxidizerY_;

            //- Fuel species
            wordList speciesFuel_;

            //- Composition of fuel [mol] mol fraction X
            map<word, scalar> fuelX_;

            //- Composition of fuel [kg] mass fraction Y
            map<word, scalar> fuelY_;

            //- Temperature of oxidizer stream
            scalar TOxidizer_{0};

            //- Temperature of fuel stream
            scalar TFuel_{0};

            //- Pressure at which the calculation take place
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


        // Debug
        const bool debug{false};


    public:

        //- Constructor
        Properties
        (
            const string&,
            const Thermo&,
            const Chemistry&
        );

        //- Destructor
        ~Properties();


        // Member functions
        
            //- Insert mfPoints_
            void insertMFPoints
            (
                const int&
            );
        
            //- Insert vmfPoints_
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
                const scalar&,
                const bool& lastEntry = false
            );
            
            //- Insert oxidizer mass composition 
            void insertCompositionOxidizerMass
            (
                const word&,
                const scalar&,
                const bool& lastEntry = false
            );
            
            //- Insert fuel mol composition 
            void insertCompositionFuelMol
            (
                const word&,
                const scalar&,
                const bool& lastEntry = false
            );
            
            //- Insert fuel mass composition 
            void insertCompositionFuelMass
            (
                const word&,
                const scalar&,
                const bool& lastEntry = false
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
            //  + word -> oxidizer (O) or fuel (F)
            void XtoY
            (
                const word 
            );

            //- Mass fraction to mol fraction
            //  + word -> oxidizer (O) or fuel (F)
            void YtoX
            (
                const word 
            );


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

            //- Return scalar dissipation rates
            scalarField sDRs() const;

            //- Return enthalpy defects
            scalarField defects() const;

            //- Return mixture fraction points
            int mfPoints() const;

            //- Return varianz of mixture fraction points
            int vmfPoints() const;

            //- Return oxidizer temperature
            scalar oxidizerTemperature() const;

            //- Return fuel temperature
            scalar fuelTemperature() const;

            //- Return runTime
            scalar runTime() const;

            //- Return deltaT
            scalar deltaT() const;

            //- Return number of defects
            unsigned int nDefects() const;

            //- Return defect value
            scalar defect
            (
                const int&
            ) const;

            //- Return pressure
            scalar p() const;

            //- Return word of input (mol or mass)
            word input() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Properties_hpp included

// ************************************************************************* //
