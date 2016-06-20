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
            
            //- Reference to thermo object
            const Thermo& thermo_;

            //- Reference to chemistry object 
            const Chemistry& chemistry_;


        // Private Data

            //- Fuel species (for adiabatic flame calculation)
            word fuel_;

            //- Oxidizer species (for adiabatic flame calculation)
            word oxidizer_;
        
            //- Mixture fraction discrete points
            int mfPoints_{0};

            //- Varianz of mixture fraction discrete points
            int vmfPoints_{0};

            //- Enthalpy defects [J/kg]
            scalarField defects_;

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
            word inertGas_{"N2"};


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

            //- ScalarDissipation rates [1/s] (stochiometric)
            scalarField sDRs_;

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
            word writeControl_{"endTime"};

            //- Write control interval
            scalar writeControlInterval_{0};

            //- DeltaT
            scalar deltaT_{0};


        // Class object


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


        // Insert functions

            //- Insert inertGas_
            void insertInertGas
            (
                const word&
            ); 

            //- Insert fuel_ 
            void insertFuel 
            (
                const word&
            ); 

            //- Insert oxidizer_ 
            void insertOxidizer
            (
                const word&
            ); 
        
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

            //- Insert interpreter keyword 
            void insertInterpreter
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
            void YtoX();

            //- Calculate enthalpy of pure fuel [J/mol]
            void HFuelAdiabatic();

            //- Calculate enthalpy of pure oxidizer [J/mol]
            void HOxidizerAdiabatic();

            //- Calculate density of pure fuel [g/m^3]
            void rhoFuel();
            
            //- Calculate density of pure oxidizer [g/m^3]
            void rhoOxidizer();

            //- Calculate defined values of pure streams
            void calcProperties();

            //- Calculate element mass fraction Zj
            map<word, scalar> calcZj
            (
                const word&  
            );

            //- Calculate Zst
            void calcZst();

            //- Calculate strain rate of stochiometric values
            void calcStrainRate();

            //- Calculate mass/mol fraction of unburned fuel and oxidizer
            //  at Zst
            void calcYXatZst();

            //- Calculate species mass fraction of burned fuel and oxidizer
            //  at Zst
            void calcYatZstBurned();

            //- Calculate adiabatic temperature at Zst [K]
            void calcAdiabaticTemperature();

            //- Calculate amount of atoms in fuel and oxidizer stream
            void calcAtomComposition();

            //- Calculate stochiometric coeffs
            void calcStochiometricCoeffs();

            //- Calculate mass and mol fraction of H2O and CO2 at Zst (burned)
            void calcBurnedSpeciesAtZst();



        // Return functions

            //- Return fuel species (for adiabatic flame calculation)
            word fuel() const;

            //- Return oxidizer species (for adiabatic flame calculation)
            word oxidizer() const;

            //- Return inert gas
            word inertGas() const;

            //- Return oxidizer species
            wordList speciesOxidizer() const;

            //- Return oxidizer mol fraction
            map<word, scalar> oxidizerX() const;
            
            //- Return oxidizer mass fraction
            map<word, scalar> oxidizerY() const;
            
            //- Return oxidizer mass fraction
            scalar oxidizerY
            (
                const word&  
            ) const;

            //- Return oxidizer element mass fraction
            map<word, scalar> oxidizerZj() const;

            //- Return fuel species
            wordList speciesFuel() const;

            //- Return fuel mol fraction
            map<word, scalar> fuelX() const;

            //- Return fuel mass fraction
            map<word, scalar> fuelY() const;

            //- Return fuel mass fraction
            scalar fuelY
            (
                const word&       
            ) const;

            //- Return fuel element mass fraction
            map<word, scalar> fuelZj() const;

            //- Return scalar dissipation rates [1/s]
            scalarField sDRs() const;

            //- Return scalar dissipation rates i [1/s]
            scalar sDRs
            (
                const int&  
            ) const;

            //- Return enthalpy defects [J/kg]
            scalarField defects() const;

            //- Return mixture fraction points
            int mfPoints() const;

            //- Return varianz of mixture fraction points
            int vmfPoints() const;

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
