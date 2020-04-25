/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     T hermo       | Open Source Thermo-Kinetic Library
  c-o-o-c     K iknetic     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Constructor.

    TKC is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or 
    (at your option) any later version.

    TKC is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TKC; if not, see <http://www.gnu.org/licenses/>

Class
    TKC::ThermoReader    

Description
    Reading the thermodynamic file and coordinate the save mechanism
    Actually, the NASA polynomials are read and set

SourceFiles
    thermoReader.cpp

\*---------------------------------------------------------------------------*/

#ifndef ThermoReader_hpp
#define ThermoReader_hpp

#include "stringManipulator.hpp"
#include "thermoData.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                      Class ThermoReader Declaration
\*---------------------------------------------------------------------------*/

class ThermoReader
:
    public StringManipulator
{
    private:

        //- List of available keywords for thermo
        wordList THERMO{"THERMO", "THERMO ALL"};

        // Private data

            // Debug
            bool debug_{false};

            //- Thermo file
            const string file_;

            //- Reference to the ThermoData object
            //  We need to access the set* functions
            ThermoData& data_;


    public:

        // Constructor and Destructor

            //- Constructor with file string and Thermo:: obj adress
            ThermoReader(const string, ThermoData&);

            //- Destructor
            ~ThermoReader();

        
        // Member functions

            //- Read the thermodynamic file
            void read();


        // Helper functions
         
            //- Find line number of keywords
            void findKeyword(int&, unsigned int&, const stringList&);

            //- Return the atomic weight of the elements 
            scalar calcWeight(const word, const scalar, const word);

            //- Construct species formula
            word constructFormula(const string);
            
            //- Split formula into species and factors
            //  + Store elements (atoms) 
            //  + Store atomic factors
            void elementsAndFactors(const string);


        // Data manipulation functions

            //- NASAPolynomial reader for first line
            void setNASAPolynomialCoeffsNo1
            (
                const string,
                const unsigned int&
            );

            //- NASAPolynomial reader for second line
            void setNASAPolynomialCoeffsNo2
            (
                const string,
                const unsigned int&
            );

            //- NASAPolynomial reader for third line
            void setNASAPolynomialCoeffsNo3
            (
                const string,
                const unsigned int&
            );

            //- NASAPolynomial reader for fourth line
            void setNASAPolynomialCoeffsNo4
            (
                const string,
                const unsigned int&
            );

            //- Calculate molecular weight of species
            void calcMolecularWeight(const word);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
