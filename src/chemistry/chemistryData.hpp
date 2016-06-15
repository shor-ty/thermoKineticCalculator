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
   AFC::ChemistryData
   
Description
    This class contains all chemistry data e.g. elements, reactions, species

SourceFiles
    chemistryData.cpp

\*---------------------------------------------------------------------------*/

#ifndef ChemistryData_hpp
#define ChemistryData_hpp

#include "typedef.hpp"
#include "math.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*--------------------------------------------------------------------------*\
                         Class ChemistryData Declaration
\*--------------------------------------------------------------------------*/

class ChemistryData
{
    private:

        // Debug switch
        bool debug{false};


        // Private data

            //- List of elements
            wordList elements_;

            //- List of species
            wordList species_;

            //- Number of dublicated reactions
            int nDuplicated_;

            //- Number of elementar reactions
            int nReac_;

            //- Reaction rate kf (forward) for each reaction
            scalarField kf_;

            //- Reaction rate kb (backward) for each reaction
            scalarList kb_;

            //- Reaction rate constant Kc for each reaction
            scalarList Kc_;

            //- Elementar reaction
            List<string> elementarReaction_;

            //- Matrix that contains all species in reaction r
            List<wordList> speciesInReaction_;

            //- Contains all reaction numbers where species i is included
            map<word, List<int>> reactionI_;

            //- boolList for TBR  
            //  + true if [M] is there
            //  + false if [M] is not included
            List<bool> TBR_;

            //- boolList for TBR Enhanced Factors
            //  + ture if [M] is modified
            //  + false if [M] represend all species
            List<bool> ENHANCE_;

            //- boolList for TBR LOW
            //  + true if 
            List<bool> LOW_;

            //- boolList for TBR TROE 
            boolList TROE_;

            //- boolList for TBR SRI 
            List<bool> SRI_;

            //- boolList for backward reaction 
            List<bool> backwardReaction_;

            //- Field that include all stochiometric values for species i
            //  in reaction r
            mapList<word, scalar> nu_;

            //- Composition of enhanced factors (species + value)
            mapList<word, scalar> enhancedFactors_;

            //- Matrix of Arrhenius coeffs
            matrix arrheniusCoeffs_;

            //- Matrix of TROE coeffs
            matrix TROECoeffs_;

            //- Matrix of Arrhenius coeffs for high pressure
            matrix LOWCoeffs_;

            //- Matrix of SRI coeffs
            matrix SRICoeffs_;

            //- Omega
            scalarField omega_;


        //- Thermodynamic available in chemistry file
        bool thermo_;


    public:

        //- Constructor
        ChemistryData();

        //- Destructor
        ~ChemistryData();


        // Member functions

            //- Set thermo_ if available in chemistry file
            void setThermo();

            //- Returnt thermo 
            bool thermo();


        // Insert functions, from ChemistryReader:: delegated 

            //- Insert elements
            void elements
            (
                const word&
            );

            //- Insert species
            void species
            (
                const word&
            );

            //- Insert elementar reaction
            void elementarReaction
            (
                const string&
            );

            //- Set arrhenius coeffs
            void arrheniusCoeffs
            (
                const scalar&,
                const scalar&,
                const scalar&
            );

            //- Insert LOW coeffs
            void LOWCoeffs
            (
                const scalar&,
                const unsigned int&
            );

            //- Insert TROE coeffs
            void TROECoeffs
            (
                const scalar&,
                const unsigned int&
            );

            //- Insert SRI coeffs
            void SRICoeffs
            (
                const scalar&,
                const unsigned int&
            );

            //- Insert ENHANCE factors (species + value) 
            void enhanceFactors
            (
                const word&,
                const scalar&
            );

            //- Increment nDuplicated_
            void incrementDuplicated();

            //- Increment nReac_
            void incrementReac();

            //- Insert species and stochiometric coeff
            void speciesNu
            (
                const word&,
                const scalar&
            );

            //- Build List<wordList> of species included in reaction r
            void speciesInReaction
            (
                const word&
            );

            //- Increment vectors and matrix size
            void incrementMatrixesVectors();


        // Setter functions, from ChemistryReader:: delegated

            //- Set the backward reaction boolean
            void BR(const bool);

            //- Set the TBR boolean
            void TBR(const bool);

            //- Set the LOW boolean
            void LOW(const bool);

            //- Set the TROE boolean
            void TROE(const bool);

            //- Set the SRI boolean
            void SRI(const bool);

            //- Set the Enhance boolean
            void ENHANCE(const bool);

            //- Set the reaction no. where species is included
            void setReacNumbers
            (
                const word&,
                const int&
            );


        // Update functions

            //- Update reaction rates kf
            void updateKf
            (
                const int&,
                const scalar&
            );
            
            //- Update reaction rates kb
            void updateKb
            (
                const int&,
                const scalar&
            );
            
            //- Update reaction rates constant Kc
            void updateKc
            (
                const int&,
                const scalar&
            );
            
            //- Update source term omega for species s
            void calculateOmega
            (
                const int&,
                const scalar&
            );
            
            //- Update source term omega for all species
            void calculateOmega
            (
                const scalarField&
            );


        // Return functions

            //- Return bool

                //- Backward reaction
                bool BR
                (
                    const int&
                ) const;

                //- Return true if elementar reaction is a TBR reaction
                bool TBR
                (
                    const int&
                ) const;

                //- Return true if elementar reaction is a LOW reaction
                bool LOW
                (
                    const int&
                ) const;

                //- Return true if elementar reaction is a TROE reaction
                bool TROE
                (
                    const int&
                ) const;

                //- Return true if elementar reaction is a SRI reaction
                bool SRI
                (
                    const int&
                ) const;

                //- Return true if elementar reaction has enhanced factors
                bool ENHANCED
                (
                    const int&
                ) const;


            //- Return all species as wordList
            wordList species() const;

            //- Return all elements as wordList
            wordList elements() const;

            //- Return no. of reaction
            int nReac() const;

            //- Return elementar reaction 
            wordList elementarReaction() const;

            //- Return elementar reaction (as string)
            word elementarReaction
            (
                const int&
            ) const;

            //- Return List of reaction no. of species
            List<int> reacNumbers
            (
                const word&
            ) const;

            //- Return species matrix for reactions
            wordMatrix speciesInReaction() const;

            //- Return species list for reaction r
            wordList speciesInReaction
            (
                const int& 
            ) const;

            //- Return arrhenius coeffs for reaction no.
            //  [0] -> pre-exponent [units depend on equation]
            //  [1] -> temperature exponent [-]
            //  [2] -> activation energy [cal/mol]
            scalarList arrheniusCoeffs
            (
                const int&
            ) const;

            //- Return arrhenius coeffs for high pressure for reaction no.
            scalarList LOWCoeffs
            (
                const int&
            ) const;

            //- Return TROE coeffs
            scalarList TROECoeffs
            (
                const int&
            ) const;

            //- Return SRI coeffs
            scalarList SRICoeffs
            (
                const int&
            ) const;

            //- Return ENHANCED species of reac no.
            wordList enhancedSpecies
            (
                const int&
            ) const;

            //- Return ENHANCED factors (species + value) of reac no.
            map<word, scalar> enhancedFactors
            (
                const int&
            ) const;

            //- Return ENHANCED factors (value) of reac no.
            scalar enhancedFactors
            (
                const int&,
                const word&
            ) const;

            //- Return reaction rate kf for reaction no.
            scalar kf
            (
                const int&
            ) const;

            //- Return reaction rates kf
            scalarList kf() const;

            //- Return reaction rate kb for reaction no.
            scalar kb
            (
                const int&
            ) const;

            //- Return reaction rates kb
            scalarList kb() const;
            
            //- Return reaction rate constant Kc for reaction no.
            scalar Kc
            (
                const int&
            ) const;

            //- Return reaction rate constant Kc
            scalarList Kc() const;
            
            //- Return omega of species s
            scalar omega
            (
                const int&
            ) const;

            //- Return omega field
            scalarField omega() const;

            //- Return map of stochiometric factors of reaction r
            map<word, scalar> nu
            (
                const int& 
            ) const;

            //- Return mapList of stochiometric factors of all reactions
            mapList<word, scalar> nu() const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ChemistryData_hpp included

// ************************************************************************* //
