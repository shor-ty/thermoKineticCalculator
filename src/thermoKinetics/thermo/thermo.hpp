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
    AFC::Thermo
    
Description
    Abstract AFC::Thermo class for thermo data and calculation

SourceFiles
    thermo.cpp

\*---------------------------------------------------------------------------*/

#ifndef Thermo_hpp
#define Thermo_hpp

#include "thermoCalc.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Thermo Declaration
\*---------------------------------------------------------------------------*/

class Thermo
:
    public ThermoCalc
{
    private:

        // Private pointer data
            
            //- Debug switch
            bool debug_{true};

            //- ThermoData object
            //ThermoData thermoData_;

            //- ThermoCalc object
            //ThermoCalc thermoCalc_;
    

    public:

        //- Constructor with fileName
        Thermo(const string, const bool thermoInChemistry = false);

        //- Destructor
        ~Thermo();


        // Insert functions delegated to thermoData

        /*
            //- Insert pressure (from properties)
            void p(const scalar);


        // Member functions
        
            //- Return species as word list
            //const wordList species() const;


            //- Moleculare weight [kg/mol]

                //- Return moleculare weight as map 
                map<word, scalar> MW() const;
            
                //- Return moleculare weight of species s
                scalar MW(const word) const;

                //- Return the mean molecular weight based on mole fraction X
                scalar MWmeanX(const map<word, scalar>&) const;


            //- Heat capacity [J/mol/K] 
            
                //- Return molar heat capacity of species s 
                //  for constant pressure
                scalar cp(const word, const scalar) const;

                //- Return heat capacity of species s 
                //  for constant volume
                scalar cv(const word, const scalar) const;


            //- Return the pressure [Pa]
            scalar p() const;


            //- Energy, enthalpy and entropy related stuff

                //- Return enthalpy of species s [J/mol]
                scalar H(const word, const scalar) const;

                //- Return entropy of species s [J/mol/K]
                scalar S(const word, const scalar) const;

                //- Return free Gibbs of species s [J/mol]
                scalar G(const word, const scalar) const;

                //- Return mean free Gibbs of species s [J/mol]
                scalar G(const scalar, const scalar, const scalar) const;

                //- Return formation enthalpy Hf == H(293) [J/mol]
                scalar Hf(const word) const;

                //- Return formation free GIBBS energy Gf = G(293) [J/mol]
                scalar Gf(const word) const;

                //- Return dHf [J/mol]
                scalar dHf(const word, const scalar) const;

                //- Return dGf [J/mol]
                scalar dGf(const word, const scalar) const;

                //- Return H0 [J/mol]
                scalar H0(const word, const scalar) const;


            //- Return the mixture concentration [mol/m^3]
            scalar C(const scalar) const;

            //- Return the density of the species [kg/m^3]
            scalar rho(const word, const scalar) const;

            //- Return low temperature bound (from NASA)
            scalar LT(const word) const;

            //- Return the phase of molecule (GAS, LIQUID, SOLID)
            word phase(const word) const;
            
            //- Return factors of atoms in species s
            scalarList elementsFactors(const word) const;


        // Summary function

            //- Build the output file that contains all data
            void summary(ostream&) const;

            //- Build NASA Coefficient table 
            void NASAPolynomials(ostream&, const word) const;

            //- Build thermoanalyse table (calcualte thermo properties)
            void thermoTable(ostream&) const;

            */
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Thermo_hpp included

// ************************************************************************* //
