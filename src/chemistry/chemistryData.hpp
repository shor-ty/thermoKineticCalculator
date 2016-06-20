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

        // Private data

            //- List of elements
            wordList elements_;

            //- List of species
            wordList species_;

            //- Number of dublicated reactions
            int nDuplicated_;

            //- Number of elementar reactions
            int nReac_;

            //- Reaction rate kf (forward)
            scalar kf_{0};

            //- Reaction rate kb (backward)
            scalar kb_{0};

            //- Reaction rate constant Kc 
            scalar Kc_{0};

            //- Thirdbody [M] concentration [g/cm^3]
            scalar M_{0};

            //- Enthalpy of reaction
            scalar dH_{0};

            //- Entropy of reaction
            scalar dS_{0};

            //- Free Gibbs energy
            scalar dG_{0};

            //- stringList for elementar reaction
            stringList elementarReaction_;

            //- Matrix that contains all species in the reaction
            wordMatrix speciesInReactions_;

            //- Species of product side of reaction r
            wordMatrix speciesInReactionProd_;

            //- Species of educt side of reaction r
            wordMatrix speciesInReactionEduc_;

            //- Contains all reaction no. where species is included
            map<word, intList> reacNumbers_;

            //- boolList for TBR  
            //  + true if [M] is there
            //  + false if [M] is not included
            boolList TBR_;

            //- boolList for TBR Enhanced Factors
            //  + ture if [M] is modified
            //  + false if [M] represend all species
            boolList ENHANCE_;

            //- boolList for TBR LOW
            //  + true if 
            boolList LOW_;

            //- boolList for TBR TROE 
            boolList TROE_;

            //- boolList for TBR SRI 
            boolList SRI_;

            //- boolList for backward reaction 
            boolList backwardReaction_;

            //- Matrix of stochiometric coeffs
            matrix nu_;

        //- Stochiometric coeffs of species of products
        vector<map<word, scalar> > nuProd_;
        
        //- Stochiometric coeffs of species of educts
        mapList<word, scalar> nuEduc_;

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
            void insertElements
            (
                const word&
            );

            //- Insert species
            void insertSpecies
            (
                const word&
            );

            //- Insert elementar reaction
            void insertElementarReaction
            (
                const string&
            );

            //- Insert arrhenius coeffs
            void insertArrheniusCoeffs
            (
                const scalar&,
                const scalar&,
                const scalar&
            );

            //- Insert LOW coeffs
            void insertLOWCoeffs
            (
                const scalar&,
                const unsigned int&
            );

            //- Insert TROE coeffs
            void insertTROECoeffs
            (
                const scalar&,
                const unsigned int&
            );

            //- Insert SRI coeffs
            void insertSRICoeffs
            (
                const scalar&,
                const unsigned int&
            );

            //- Insert ENHANCE factors (species + value) 
            void insertEnhanceFactors
            (
                const word&,
                const scalar&
            );

            //- Increment nDuplicated_
            void incrementDuplicated();

            //- Increment nReac_
            void incrementReac();

            //- Insert stochiometric coeffs
            void insertNu
            (
                const scalar&
            );

            //- Insert stochiometric coeffs for products and species
            void insertProd
            (
                const word&,
                const scalar&
            );

            //- Insert stochiometric coeffs for educts and species
            void insertEduc
            (
                const word&,
                const scalar&
            );

            //- Insert species (reac|prod) 
            void insertReacProd
            (
                const word&
            );

            //- Increment vectors and matrix size
            void incrementMatrixesVectors();


        // Setter functions, from ChemistryReader:: delegated

            //- Set the backward reaction boolean
            void setBR();

            //- Set the TBR boolean
            void setTBR();

            //- Set the LOW boolean
            void setLOW();

            //- Set the TROE boolean
            void setTROE();

            //- Set the SRI boolean
            void setSRI();

            //- Set the Enhance boolean
            void setENHANCE();

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
                const scalar&
            );
            
            //- Update reaction rates kb
            void updateKb
            (
                const scalar&
            );
            
            //- Update reaction rates constant Kc
            void updateKc
            (
                const scalar&
            );

            //- Update dH
            void updateDH
            (
                const scalar&
            );

            //- Update dS
            void updateDS
            (
                const scalar&
            );

            //- Update dG
            void updateDG
            (
                const scalar&
            );
            
            //- Update M
            void updateM
            (
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

                //- LOW
                bool LOW
                (
                    const int&
                ) const;

                //- TROE
                bool TROE
                (
                    const int&
                ) const;

                //- TBR
                bool TBR
                (
                    const int&
                ) const;

                //- SRI
                bool SRI
                (
                    const int&
                ) const;

                //- ENHANCED
                bool ENHANCED
                (
                    const int&
                ) const;

            //- Return M_ [g/cm^3]
            scalar M() const;

            //- Return dH
            scalar dH() const;

            //- Return dH
            scalar dS() const;

            //- Return dG
            scalar dG() const;

            //- Return all elements
            wordList elements() const;

            //- Return all species
            wordList species() const;

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
            intList reacNumbers
            (
                const word&
            ) const;

            //- Return species matrix for reactions
            wordMatrix speciesInReactions() const;

            //- Return species list for reaction r
            wordList speciesInReaction
            (
                const int& 
            ) const;

            //- Return species list that act as product in reaction r
            wordList prodSpecies
            (
                const int&
            ) const;
            
            //- Return species list that act as educt in reaction r
            wordList educSpecies
            (
                const int&
            ) const;

            //- Return arrhenius coeffs for reaction no.
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

            //- Return reaction rate kf 
            scalar kf() const;

            //- Return reaction rate kb 
            scalar kb() const;

            //- Return reaction rate constant Kc 
            scalar Kc() const;
            
            //- Return omega of species s
            scalar omega
            (
                const int&
            ) const;

            //- Return omega field
            scalarField omega() const;

            //- Return stochiometric factors of reaction r
            scalarList nu
            (
                const int&
            ) const;

            //- Return stochiometric factors of products of reaction r
            map<word, scalar> nuProd
            (
                const int&
            ) const;

            //- Return stochiometric factors of products of reaction r
            map<word, scalar> nuEduc
            (
                const int&
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ChemistryData_hpp included

// ************************************************************************* //
