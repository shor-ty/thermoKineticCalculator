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

            //- stringList for elementar reaction
            stringList elementarReaction_;

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
            boolList kb_;

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

            //- Increment vectors and matrix size
            void incrementMatrixesVectors();


        // Setter functions, from ChemistryReader:: delegated

            //- Set the backward reaction boolean
            void setKB();

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

        
        // Return functions

            //- Return all species as wordList
            wordList species() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ChemistryData_hpp included

// ************************************************************************* //
