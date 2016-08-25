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

        // Private class data

        // Debug
        const bool debug_{false};

        // Private Data

            //- Fuel species (for adiabatic flame calculation)
            word fuel_;

            //- Oxidizer species (for adiabatic flame calculation)
            word oxidizer_;
        
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

            //- Composition of oxidizer element mass fraction Zj [-]
            map<word, scalar> oxidizerZj_;

            //- Atoms and amount in oxidizer
            map<word, scalar> oxidizerA_;

            //- Fuel species
            wordList speciesFuel_;

            //- Composition of fuel mol fraction X [-]
            map<word, scalar> fuelX_;

            //- Composition of fuel mass fraction Y [-]
            map<word, scalar> fuelY_;

            //- Composition of fuel element mass fraction Zj [-]
            map<word, scalar> fuelZj_;

            //- Atoms and amount in fuel
            map<word, scalar> fuelA_;

            //- Temperature of oxidizer stream [K]
            scalar TOxidizer_{0};

            //- Temperature of fuel stream [K]
            scalar TFuel_{0};

            //- Pressure at which the calculation take place [Pa]
            scalar p_{0};

            //- For interpreter
            word interpreter_{"FALSE"};

            //- Inert gas
            word fuelSpecies_;

            //- Inert gas
            word oxidizerSpecies_;

            //- Inert gas
            word inertSpecies_{"N2"};


        // Data at stochiometric condition 

            //- Composition of mass fraction Y at Zst (unburned) [-]
            map<word, scalar> YatZstu_;

            //- Composition of mol fraction X at Zst (unburned) [-]
            map<word, scalar> XatZstu_;

            //- Composition of mass fraction Y at Zst (burned) [-]
            map<word, scalar> YatZstb_;

            //- Composition of mol fraction X at Zst (burned) [-]
            map<word, scalar> XatZstb_;

            //- Composition of element mass fraction Zj at Zst [-]
            //  unburned == burned
            map<word, scalar> ZjatZst_;

            //- Adiabatic flame temperature [K] (simplified)
            scalar Tadiabatic_{0};

            //- Stochiometric mixture fraction Zst [-]
            scalar Zst_{0};

            //- Strain rates [1/s] 
            scalarField as_;

            //- Stochiometric coeff for CO2
            scalar nuCO2_{0};

            //- Stochiometric coeffs for H2O
            scalar nuH2O_{0};

            //- Stochiometric coeffs for O2
            scalar nuO2_{0};

            //- Min O2 for combustion
            scalar omin_{0};


        // Constant data
           
            //- Adiabatic enthalpy of pure fuel [J/kg]
            scalar fuelH_{0};
           
            //- Adiabatic enthalpy of pure oxidizer [J/kg]
            scalar oxidizerH_{0};

            //- Density of pure fuel 
            scalar fuelRho_{0};

            //- Density of pure oxidizer
            scalar oxidizerRho_{0};


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
            //      + time
            word writeControl_{"time"};

            //- Write control interval
            scalar writeControlInterval_{0};

            //- Write control time
            scalar writeControlTime_{0};

            //- Deltat
            scalar deltat_{0};

            //- Current time 
            scalar currentTime_{0};


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
            //  TODO sort after everything is read
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
            
            //- Insert write control time (save each x seconds)
            void insertWriteControlTime
            (
                const scalar&
            );

            //- Insert deltat
            void insertDeltat
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

            //- Insert interpreter keyword 
            void insertInterpreter
            (
                const word& 
            );

            //- Insert fuel species
            void fuelSpecies
            (
                const word&
            );

            //- Set oxidizer species
            void oxidizerSpecies
            (
                const word&
            );

            //- Set inert species
            void inertSpecies
            (
                const word&
            );


        // Other functions

            //- Initial boundarys with all species = 0
            void initialBoundary();

            //- Check if all data are set
            void check();

            //- Convert mol or mass fraction of pure streams
            void convertFractions();

            //- Mol fraction to mass fraction
            void XtoY();

            //- Mass fraction to mol fraction
            //  + word -> oxidizer (O) or fuel (F)
            /*void YtoX
            (
                const word 
            );*/

        // Time RelateD Functions 

            //- Update the current time stamp [s]
            void updateCurrentTime
            (
                const scalar&
            );

            //- Return runTime [s]
            scalar runTime() const;

            //- Return deltat
            scalar deltat() const;

            //- Return the current time [s]
            scalar currentTime() const;

            //- Return write control for time
            scalar write() const;


        // Return Functions

            //- Return fuel species (for adiabatic flame calculation)
            word fuel() const;

            //- Return oxidizer species (for adiabatic flame calculation)
            word oxidizer() const;

            //- Return inert species
            word inert() const;

            //- Return species of oxidizer (Z = 0)
            wordList speciesOxidizer() const;

            //- Return mol fraction of oxidizer (Z = 0)
            map<word, scalar> oxidizerX() const;

            //- Return mass fraction of oxidizer (Z = 0)
            map<word, scalar> oxidizerY() const;

            //- Return element mass fraction of oxidizer (Z=0) 
            map<word, scalar> oxidizerZj() const;

            //- Return mof fraction of species s in oxidizer (Z = 0)
            scalar oxidizerX
            (
                const word&
            ) const;

            //- Return mass fraction of species s in oxidizer (Z = 0)
            scalar oxidizerY
            (
                const word&
            ) const;

            //- Return species of fuel (Z = 1)
            wordList speciesFuel() const;

            //- Return mol fraction of fuel (Z = 1)
            map<word, scalar> fuelX() const;

            //- Return mass fraction of fuel (Z = 1)
            map<word, scalar> fuelY() const;

            //- Return element mass fraction of fuel (Z = 1)
            map<word, scalar> fuelZj() const;

            //- Return mol fraction of species s of fuel (Z = 1)
            scalar fuelX
            (
                const word&
            ) const;

            //- Return mass fraction of species s of fuel (Z = 1)
            scalar fuelY
            (
                const word&
            ) const;

            //- Return scalar dissipation rates [1/s]
            scalarField sDRs() const;

            //- Return scalar dissipation rate i [1/s]
            scalar sDRs
            (
                const int&  
            ) const;

            //- Return enthalpy defects [J/kg]
            scalarField defects() const;

            //- Return the number of discrete points of Z
            int nZPoints() const;

            //- Return the number of discrete points of Zvar
            int nZvarPoints() const;

            //- Return the oxidizer temperature [K] (Z = 0)
            scalar oxidizerTemperature() const;

            //- Return the fuel temperature [K] (Z = 1)
            scalar fuelTemperature() const;

            //- Return the number of defects 
            unsigned int nDefects() const;

            //- Return the enthalpy defect value [J/kg]
            scalar defect
            (
                const int&
            ) const;

            //- Return pressure [Pa]
            scalar p() const;

            //- Return word of input (mol or mass)
            word input() const;

            //- Return if thermo data or chemistry data are analysed
            word interpreter() const;

            //- Return Zst
            scalar Zst() const;

            //- Return mass fraction at stochiometric mixture fraction Zst [-]
            //  unburned state
            scalar YatZstu
            (
                const word&
            ) const;

            //- Return mass fraction at stochiometric mixture fraction Zst [-]
            //  burned state
            scalar YatZstb
            (
                const word&
            ) const;

            //- Return adiabatic flame temperature [K]
            scalar Tadiabatic() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Properties_hpp included

// ************************************************************************* //
