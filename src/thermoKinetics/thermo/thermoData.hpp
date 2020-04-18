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
   AFC::ThermoData
   
Description
    This class contains all thermo data e.g. coeffs of NASA polynomials,
    the pressure, species names (nick-names and chemical formula),...

SourceFiles
    thermoData.cpp

\*---------------------------------------------------------------------------*/

#ifndef ThermoData_hpp
#define ThermoData_hpp

#include "definitions.hpp"

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

            //- Debug switch
            bool debug_{false};

            //- Thermodynamic is located in the chemistry file
            bool thermoInChemistry_;

            
            //- Pressure [Pa]
            scalar p_;


            //- Data related to species
            
                //- List of species names -> ACETOL
                wordList species_;

                //- List of species formula -> C3H6O2
                wordList formula_;

                //- List of elements in each species
                map<word, wordList> elementsInSpecies_;

                //- List of stochiometric factors of each element in species
                map<word, scalarList> elementAtoms_;

                //- Elements and amount in species (Summary of both above)
                map<word, map<word, scalar> > elements_;

                //- Molecular weight of each single species [kg/mol]
                map<word, scalar> MW_;

                //- Phase status of each species
                map<word, word> phase_;


            //- NASA Polynomial related information

                //- List of Low Temperature limit for NASA of each species
                map<word, scalar> LT_;

                //- List of Common Temperature limit for NASA of each species
                map<word, scalar> CT_;

                //- List of High Temperature limit for NASA of each species
                map<word, scalar> HT_;

                //- List of the polynomial coefficients for LT of each species
                map<word, scalarField> NASACoeffsLT_;

                //- List of the polynomial coefficients for HT of each species
                map<word, scalarField> NASACoeffsHT_;


    public:


        //- Constructor that take the path to the file one wants to read
        //  The second argument is related to cases where the thermodynamic
        //  is located inside the chemistry (kinetic) file        
        ThermoData(const string, const bool thermoInChemistry = false);

        //- Destructor
        ~ThermoData();


        // Insert functions, from Thermo:: delegated

            //- Insert the pressure
            void p(const scalar);


        // Insert functions, from ThermoReader:: delegated 

            //- Insert species
            void setSpecies(const word);

            //- Insert chemical species (formula)
            void setChemicalFormula(const word);

            //- Insert element and factor of species
            void setElementAndAtoms(const word, const unsigned int);

            //- Insert molecular weight [kg/mol]
            void setMolecularWeight(const scalar);

            //- Insert phase of species
            void setPhase(const word);

            //- Insert low temperature
            void setLT(const scalar);

            //- Insert high temperature
            void setHT(const scalar);

            //- Insert common temperature
            void setCT(const scalar);

            //- Insert hight temperature poly coeffs
            void setNASACoeffsHT(const scalar);

            //- Insert low temperature poly coeffs
            void setNASACoeffsLT(const scalar);

            //- Update the elements map (stores elements and factors)
            void updateElementsAndFactors();


        // Return functions

            //- Return the pressure [Pa]
            scalar p() const;

            //- Return species as wordList
            const wordList species() const;

            //- Return the formula of species as wordList
            const wordList formula() const;

            //- Return the elements of species s
            const wordList elementsInSpecies(const word) const;
            
            //- Return the elements of species (chemical form) s 
            const wordList elementsInSpeciesChem(const word) const;

            //- Return the factor of elements in species s
            const scalarList elementAtoms(const word) const;

            //- Return the factor of elements as map of species s
            const map<word, scalar> elementAtomsMap(const word) const;

            //- Return the factor of elements in species (chemical form) s 
            const map<word, scalar> elementAtomsChem(const word) const;
            
            //- Return moleculare weight as map [g/mol]
            const map<word, scalar> MW() const;

            //- Return moleculare weight of species s [g/mol]
            scalar MW(const word) const;

            //- Return the phase of species as map
            const map<word, word> phase() const;

            //- Return the phase of species s
            const word phase(const word) const;
            
            //- Return LOW temperature of polynomials of species s
            scalar LT(const word) const;

            //- Return COMMON temperature of polynomials of species s
            scalar CT(const word) const;

            //- Return HIGH temperature of polynomials of species s
            scalar HT(const word) const;

            //- Return polyCoeffs for HIGH temperature
            const List<scalar> NASACoeffsHT(const word) const;

            //- Return polyCoeffs for LOW temperature
            const List<scalar> NASACoeffsLT(const word) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ThermoData_hpp included

// ************************************************************************* //
