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
   AFC::ThermoData
   
Description
    This class contains all thermo data e.g. elements, reactions, species

SourceFiles
    thermoData.cpp

\*---------------------------------------------------------------------------*/

#ifndef ThermoData_hpp
#define ThermoData_hpp

#include "typedef.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*--------------------------------------------------------------------------*\
                         Class ThermoData Declaration
\*--------------------------------------------------------------------------*/

class ThermoData
{
    private:

        // Private data

            //- Pressure for calculation
            scalar p_;

            //- List of species (official species name)
            wordList species_;

            //- List of species formula (chemical formula)
            //  + ACETOL -> C3H6O2
            wordList formula_;

            //- map of List of atoms of each species
            // wordMatrix speciesAtoms_;
            map<word, wordList> speciesAtoms_;

            //- map of List of coeffs of atoms
            // matrix atomFactors_;
            map<word, scalarList> atomFactors_;

            //- Hashtable of molecular weight of species
            map<word, scalar> MW_;

            //- Hashtable of phase of species
            map<word, word> phase_;

            //- Hashtable of low temperature
            map<word, scalar> LT_;

            //- Hashtable of high temperature
            map<word, scalar> HT_;

            //- Hashtable of common temperature
            map<word, scalar> CT_;

            //- Hashtable of polycoeffs for high temperature range
            map<word, scalarField> NASACoeffsHT_;

            //- Hashtable of polycoeffs for high temperature range
            map<word, scalarField> NASACoeffsLT_;


        //- Thermodynamic available in chemistry file
        bool thermo_;


    public:

        //- Constructor that take ChemistryData::thermo_
        ThermoData
        (
            const bool& 
        );

        //- Destructor
        ~ThermoData();


        // Member functions



        // Insert functions, from ThermoReader:: delegated 

            //- Insert species
            void insertSpecies
            (
                const word&
            );

            //- Insert chemical species (formula)
            void insertChemicalFormula
            (
                const word&
            );

            //- Insert atom and factor of species
            void insertAtomAndFactor
            (
                const word&,
                const unsigned int&
            );

            //- Insert molecular weight
            void insertMolecularWeight
            (
                const scalar&
            );

            //- Insert phase of species
            void insertPhase
            (
                const word&
            );

            //- Insert low temperature
            void insertLT
            (
                const scalar&
            );

            //- Insert high temperature
            void insertHT
            (
                const scalar&
            );

            //- Insert common temperature
            void insertCT
            (
                const scalar&
            );

            //- Insert hight temperature poly coeffs
            void insertNASACoeffsHT
            (
                const scalar&
            );

            //- Insert low temperature poly coeffs
            void insertNASACoeffsLT
            (
                const scalar&
            );


        // Insert functions, from Thermo:: delegated

            //- Insert the pressure
            void p
            (
                const scalar&
            );


        // Setter functions, from ThermoReader:: delegated

        // Return functions

            //- Return species as wordList
            wordList species() const;

            //- Return moleculare weight of species s
            scalar MW
            (
                const word&
            ) const;

            //- Return moleculare weight as map
            map<word, scalar> MW() const;

            //- Return LOW temperature of polynomials of species s
            scalar LT
            (
                const word&
            ) const;

            //- Return COMMON temperature of polynomials of species s
            scalar CT 
            (
                const word&
            ) const;

            //- Return HIGH temperature of polynomials of species s
            scalar HT 
            (
                const word&
            ) const;

            //- Return polyCoeffs for HIGH temperature
            scalarField NASACoeffsHT 
            (
                const word&
            ) const;

            //- Return polyCoeffs for LOW temperature
            scalarField NASACoeffsLT 
            (
                const word&
            ) const;

            //- Return the pressure [Pa]
            scalar p() const;

            //- Return the atoms of species
            wordList speciesAtoms
            (
                const word&
            ) const;
            
            //- Return the atoms of species (chemical form)
            //  Not implemented
            wordList speciesAtomsChem
            (
                const word&
            ) const;

            //- Return the factor of atoms in species
            scalarList atomFactors
            (
                const word&
            ) const;

            //- Return the factor of atoms as map
            map<word, scalar> atomFactors() const;

            //- Return the factor of atoms in species (chemical form)
            //  Not implemented
            scalarList atomFactorsChem
            (
                const word&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ThermoData_hpp included

// ************************************************************************* //
