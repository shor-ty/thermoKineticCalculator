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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

// Forward declaration
class ChemistryReader;

/*---------------------------------------------------------------------------*\
                            Class Chemistry Declaration
\*---------------------------------------------------------------------------*/

class Chemistry
{
    private:

        // Private pointer data

            //- Pointer to ChemistryReader object 
            smartPtr<ChemistryReader> pCR_;

            //- Pointer to ChemistryData object
            smartPtr<ChemistryData> pCD_;
    

    public:

        //- Constructor
        Chemistry();

        //- Destructor
        ~Chemistry();


        // Runtime object creator

            //- Generate new objects
            void New
            (
                const string&
            );


        // Member Functions
            
            //- read chemistry file
            void readChemistry();


        // Insert functions from ChemistryReader:: (aggregation)
        // -> delegate to ChemitryData::

            //- Insert elements -> delegated
            void insertElements
            (
                const word&
            );

            //- Insert species -> delegated
            void insertSpecies
            (
                const word&
            );

            //- Insert elementar reaction -> delegated
            void insertElementarReaction
            (
                const string& 
            );

            //- Insert arrhenius coeffs -> delegated
            void insertArrheniusCoeffs
            (
                const scalar&,
                const scalar&,
                const scalar&
            );

            //- Insert LOW coeffs -> delegated
            void insertLOWCoeffs
            (
                const scalar&,
                const unsigned int&
            );

            //- Insert TROE coeffs -> delegated
            void insertTROECoeffs
            (
                const scalar&,
                const unsigned int&
            );
            
            //- Insert SRI coeffs -> delegated
            void insertSRICoeffs
            (
                const scalar&,
                const unsigned int&
            );

            //- Insert ENHANCE value -> delegate
            void insertMvalue
            (
                const scalar&
            );

            //- Insert ENHANCE comp -> delegate
            void insertMcomp
            (
                const string&
            );

            //- Increment nDuplicated_ -> delegated
            void incrementDuplicated();

            //- Increment nReac_ -> delegated
            void incrementReac();
                
            //- Increment vectors and matrix size -> delegated
            void incrementMatrixesVectors();


        // Setter bool functions

            //- Set the backward reaction boolean -> delegated
            void setKB();

            //- Set the TBR boolean -> delegated
            void setTBR();

            //- Set the LOW boolean -> delegated
            void setLOW();

            //- Set the TROE boolean -> delegated
            void setTROE();

            //- Set the SRI boolean -> delegated
            void setSRI();

            //- Set the Enhance boolean --> delegated
            void setENHANCE(); 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Chemistry_hpp included

// ************************************************************************* //
