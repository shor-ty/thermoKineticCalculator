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

            //- Reaction rate kf (forward) for each reaction
            scalarField kf_;

            //- Reaction rate kb (backward) for each reaction
            scalarList kb_;

            //- Reaction rate constant Kc for each reaction
            scalarList Kc_;

            //- stringList for elementar reaction
            stringList elementarReaction_;

            //- Matrix that contains all reactions and for each
            //  reaction all included species
            wordMatrix speciesInReactions_;

            //- Contains all reaction no. where species i is included
            matrix reacNoSpecies_;

            //- boolList for TBR  
            boolList TBR_;

            //- boolList for TBR LOW
            boolList LOW_;

            //- boolList for TBR TROE 
            boolList TROE_;

            //- boolList for TBR SRI 
            boolList SRI_;

            //- boolList for TBR Enhanced Factors
            boolList ENHANCE_;

            //- boolList for backward reaction 
            boolList backwardReaction_;

            //- Matrix of stochiometric coeffs
            matrix nu_;

            //- Matrix of TB M (composition of species)
            wordMatrix Mcomp_;

            //- Matrix of TB M (values of species)
            matrix Mvalue_;

            //- Matrix of Arrhenius coeffs
            matrix arrheniusCoeffs_;

            //- Matrix of TROE coeffs
            matrix TROECoeffs_;

            //- Matrix of Arrhenius coeffs for LOW pressure
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

            //- Insert ENHANCE value
            void insertMvalue
            (
                const scalar&
            );

            //- Insert ENHANCE comp
            void insertMcomp
            (
                const string&
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

            //- Set the reaction no. (for species)
            void setReacNoSpecies
            (
                const int&,
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
            void updateOmega
            (
                const int&,
                const scalar&
            );
            
            //- Update source term omega for all species
            void updateOmega
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


            //- Return all species as wordList
            wordList species() const;

            //- Return no. of reaction
            int nReac() const;

            //- Return elementar reaction (as string)
            word elementarReaction
            (
                const int&
            ) const;

            //- Return scalarList of reaction nos. for species
            scalarList reacNoForSpecies
            (
                const int&
            ) const;

            //- Return species matrix for reactions
            wordMatrix speciesInReactions() const;

            //- Return species list for reaction r
            wordList speciesInReaction
            (
                const int& 
            ) const;

            //- Return arrhenius coeffs for reaction no.
            scalarList arrheniusCoeffs
            (
                const int&
            ) const;

            //- Return Mcomp for reaction no.
            wordList Mcomp
            (
                const int&
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

            //- Return stochiometric factors of reaction r
            scalarList nu
            (
                const int&
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ChemistryData_hpp included

// ************************************************************************* //
