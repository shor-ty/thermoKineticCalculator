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
    AFC::Chemistry
    
Description
    Abstract AFC::Chemistry class for chemistry data and calculation

SourceFiles
    chemistry.cpp

\*---------------------------------------------------------------------------*/

#ifndef Chemistry_hpp
#define Chemistry_hpp

#include "chemistryReader.hpp"
#include "chemistryData.hpp"
#include "chemistryCalc.hpp"
#include "thermo.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Chemistry Declaration
\*---------------------------------------------------------------------------*/

class Chemistry
{
    public:

        // Private reference data

            //- ChemistryData obj
            ChemistryData chemData_;

            //- ChemistryCalc obj
            ChemistryCalc chemCalc_;


        // Debug 
        bool debug{false};

    public:

        //- Constructor
        Chemistry
        (
            const string&
        );

        //- Destructor
        ~Chemistry();


        // Member Functions
            
            //- Read chemistry file
            void readChemistry();

            //- Return chemistryData::thermo_
            bool thermo();


        // Calculation Functions

            //- Calculate the source term of each species (omega)
            scalar calculateOmega
            (
	            const word&,
                const scalar&,
                map<word, scalar>&,
                const Thermo&
            );

            //- Calculate forward reaction rate kf ::Interpreter
            void calculateKf
            (
                const int&,
                const scalar&  
            );

            //- Calculate equilibrium reaction rate Kc ::Interpreter
            void calculateKc
            (
                const int&,
                const scalar&,
                const Thermo&
            );

            //- Calculate backward reaction kb ::Interpreter
            void calculateKb ();


        // Update Functions

            //- Update [M] (only inert gas) ::Interpreter
            void updateM
            (
                const scalar&
            );

        // Create Functions

            //- Create field that contains all reaction no. for each species
            void createSpeciesInReaction();

            //- Insert reaction no. for species ::DELEGATED to CHEMISTRYDATA
            void insertReacNoForSpecies
            (
                const int&,
                const int&
            );


        // Return Functions

            //- Return elements
            wordList elements() const;
        
            //- Retrun species
            wordList species() const;

            //- Return kf
            scalar kf() const;

            //- Return Kc
            scalar Kc() const;

            //- Return kb
            scalar kb() const;

            //- Return no. of elementar reaction
            //int nReac() const;

            //- Return elementar reaction (as string) DELEGATED
            word elementarReaction
            (
                const int&
            ) const;

            //- Return all elementar reactions
            stringList elementarReaction() const;

            //- Return arrhenius coeffs ::Interpreter
            scalarList arrheniusCoeffs
            (
                const int&
            ) const;

            //- Return LOW arrhenius coeffs ::Interpreter
            scalarList LOWCoeffs
            (
                const int&
            ) const;

            //- Return TROE coeffs ::Interpreter
            scalarList TROECoeffs
            (
                const int&
            ) const;

            //- Return SRI coeffs ::Interpreter
            scalarList SRICoeffs
            (
                const int&
            ) const;

            //- Return TBR bool ::Interpreter
            bool TBR
            (
                const int&
            ) const;

            //- Return ENHANCED bool ::Interpreter
            bool ENHANCED
            (
                const int&
            ) const;

            //- Return LOW bool ::Interpreter
            bool LOW
            (
                const int&
            ) const;

            //- Return TROE bool ::Interpreter
            bool TROE
            (
                const int&
            ) const;

            //- Return SRI bool ::Interpreter
            bool SRI 
            (
                const int&
            ) const;

            //- Return enhanced factors ::Interpreter
            map<word, scalar> enhancedFactors
            (
                const int&
            ) const;

            //- Return enhanced species ::Interpreter
            wordList enhancedSpecies
            (
                const int&
            ) const;

            //- Return enthalpy of actual reaction
            scalar dH() const;

            //- Return entropy of actual reaction
            scalar dS() const;

            //- Return free Gibbs energy of actual reaction
            scalar dG() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Chemistry_hpp included

// ************************************************************************* //
