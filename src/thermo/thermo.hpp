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
    AFC::Thermo
    
Description
    Abstract AFC::Thermo class for thermo data and calculation

SourceFiles
    thermo.cpp

\*---------------------------------------------------------------------------*/

#ifndef Thermo_hpp
#define Thermo_hpp

#include "thermoReader.hpp"
#include "thermoData.hpp"
#include "thermoCalc.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Thermo Declaration
\*---------------------------------------------------------------------------*/

class Thermo
{
    private:

        // Private pointer data

            //- ThermoData object
            ThermoData thermoData_;

            //- ThermoCalc object
            ThermoCalc thermoCalc_;
    

    public:

        //- Constructor with fileName
        Thermo
        (
            const string&,
            const bool&
        );

        //- Destructor
        ~Thermo();


        // Insert functions delegated to thermoData

            //- Insert pressure (from properties)
            void p
            (
               const scalar&
            );


        // Return | calculation functions
        
            //- Return species as word list
            wordList species() const;

            //- Return moleculare weight of species s [g/mol]
            scalar MW
            (
                const word&
            ) const;

            //- Return moleculare weight as map 
            map<word, scalar> MW() const;

            //- Return the mean molecular weight based on X [g/mol]
            scalar MmeanX
            (
                const map<word, scalar>&
            ) const;

            //- Return heat capacity of species s [J/mol/K]
            scalar cp
            (
                const word&,
                const scalar&
            ) const;

            //- Return enthalpy of species s [J/mol]
            scalar H
            (
                const word&,
                const scalar&
            ) const;

            //- Return entropy of species s [J/mol/K]
            scalar S
            (
                const word&,
                const scalar&
            ) const;

            //- Return free Gibbs of species s [J/mol]
            scalar G
            (
                const word&,
                const scalar&
            ) const;

            //- Return mean free Gibbs of species s [J/mol]
            scalar G
            (
                const scalar&,
                const scalar&,
                const scalar&
            ) const;

            //- Return the pressure [Pa]
            scalar p() const;

            //- Return formation enthalpy Hf == H(293) [J/mol]
            scalar Hf
            (
                const word&
            ) const;

            //- Return formation free GIBBS energy Gf = G(293) [J/mol]
            scalar Gf
            (
                const word&
            ) const;

            //- Return dH [J/mol]
            scalar dH
            (
                const word&,
                const scalar&
            ) const;

            //- Return dG [J/mol]
            scalar dG
            (
                const word&,
                const scalar&
            ) const;

            //- Return the mixture concentration [mol/m^3]
            scalar C
            (
                const scalar&
            ) const;

            //- Return low temperature bound (from NASA)
            scalar LT
            (
                const word&
            ) const;

            //- Return the phase of molecule (GAS, LIQUID, SOLID)
            word phase
            (
                const word&
            ) const;

        // Summary function

            void summary
            (
                ostream&
            ) const;

            //- Return factors of atoms in species s
            scalarList elementsFactors
            (
                const word&
            ) const;

            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Thermo_hpp included

// ************************************************************************* //
