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
    AFC::ChemistryReader

Description
    Reading the chemkin III file

SourceFiles
    chemistryReader.cpp

\*---------------------------------------------------------------------------*/

#ifndef ChemistryReader_hpp
#define ChemistryReader_hpp

#include "stringManipulator.hpp"
#include "chemistryData.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                      Class ChemistryReader Declaration
\*---------------------------------------------------------------------------*/

class ChemistryReader
:
    public StringManipulator
{
    private:


        //- List of available keywords for element block
        wordList ELEMENT{"ELEMENTS", "ELEM"};

        //- List of available keywords for species block
        wordList SPECIES{"SPECIES", "SPEC"};

        //- List of available keywords for thermo block
        wordList THERMO{"THERMO", "THERMO ALL"};

        //- List of available keywords for reaction block
        wordList REACTION{"REACTIONS", "REACT"};


        // Private data

            // Debug switch
            bool debug_{false};

            //- Chemistry file
            string file_;


    public:

        // Constructor and Destructor

            //- Constructor with file string and Chemistry:: obj adress
            ChemistryReader(const string);

            //- Destructor
            ~ChemistryReader();


        // Member functions

            //- Read chemistry file and delegate data
            void read(ChemistryData&);

            //- Take all data from ELEMENT block
            void readElementBlock(const stringList&, ChemistryData&);

            //- Take all data from SPECIES block
            void readSpeciesBlock(const stringList&, ChemistryData&);

            //- Take all data from THERMO block
            void readThermoBlock(const stringList&, ChemistryData&);

            //- Take all data from REACTION block
            void readReactionBlock(const stringList&, ChemistryData&);


        // Helper functions

            //- Find line number of keyword
            void findKeyword
            (
                int&,
                unsigned int&,
                const stringList&,
                const string
            );

            //- Return string between '/' and '/'
            stringList extractData(const string);


        // Data manipulation functions

            //- Manipulate elementar reaction string and analyze reaction
            void analyzeReaction(const stringList&, ChemistryData&);

            //- Manipulate LOW coeffs
            void LOWCoeffs(const string, const unsigned int, ChemistryData&);

            //- Manipulate TROE coeffs
            void TROECoeffs(const string, const unsigned int, ChemistryData&);

            //- Manipulate SRI coeffs
            void SRICoeffs(const string, const unsigned int, ChemistryData&);

            //- Manipulate ENHANCE factors
            void enhanceFactors(const string, ChemistryData&);

            //- Analyse reaction site (product or reactants)
            void analyzeReacSite(string, const word, ChemistryData&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
