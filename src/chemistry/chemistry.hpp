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
    private:

        // Private reference data

            //- ChemistryData obj
            ChemistryData chemData_;

            //- ChemistryCalc obj
            const ChemistryCalc chemCalc_;

            //- Thermo obj
            const Thermo& thermo_;


        // Debug 
        bool debug_{false};

    public:

        //- Constructor
        Chemistry
        (
            const string&,
            const Thermo&
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

            //- Calculate forward reaction rate kf
            scalar kf
            (
                const int&,
                const scalar&,
                const bool LOW = false
            ) const;

            //- Calculate backward reaction rate kf
            scalar kb
            (
                const int&,
                const scalar&,  
                const bool LOW = false
            ) const;

            //- Calculate equilibrium reaction rate keq
            scalar keq 
            (
                const int&,
                const scalar&
            ) const;

            //- Calculate Fcent for TROE formula
            scalar Fcent
            (
                const int&,
                const scalar&
            ) const;

            //- Calculate logF for TROE formula
            scalar Flog
            (
                const int&,
                const scalar&,
                const scalar&
            ) const;


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

            //- Insert reaction no. for species delegated to CHEMISTRYDATA
            void insertReacNoForSpecies
            (
                const int&,
                const int&
            );

        // Return Functions

            //- Return bool

                //- Return true if elementar reaction is handled with
                //  forward and backward reaction character
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



            //- Retrun all species that are included in the reactions
            wordList species() const;

            //- Return all elements that are included in the reactions
            wordList elements() const;

            unsigned int nDublicated() const;

            //- Return no. of elementar reaction
            int nReac() const;

            //- Return elementar reaction r (as string) 
            string elementarReaction
            (
                const int&
            ) const;

            //- Return all elementar reaction (as wordList)
            List<string> elementarReaction() const;

            //- Return all reaction no of species 
            /*scalarField reacNoForSpecies
            (
                const int&
            ) const;

            //- Return reaction rates k
            scalarField k() const;

            //- Return wordMatrix for species in reactions
            wordMatrix speciesInReactions() const;

            //- Return reaction rate k of reaction no.
            scalar k
            (
                const int& 
            ) const;*/


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

            //- Return enthalpy of reaction r and temperature T
            scalar dH
            (
                const int&,
                const scalar&
            ) const;
            
            //- Return free GIBBS energy of reaction r and temperature T
            scalar dG
            (
                const int&,
                const scalar&
            ) const;

            //- Return entropy of reaction r and temperature T
            scalar dS
            (
                const int&,
                const scalar&
            ) const;


        //- Summary Functions

            //- Build the summary as ostream
            void summary
            (
                ostream&   
            ) const;

            //- Build chemical table for reaction r (for summary)
            void chemicalTable
            (
                ostream&
            ) const;

            //- Build table for kf and kb based on general coefficients
            //  or if we specify "LOW" then based on low coefficients
            void buildTablekf
            (
                const int&,
                ostream&,
                const bool LOW = false 
            ) const;

            //- Build TROE table (calculate F_cent and logF)
            void buildTROETable
            (
                const int&,
                ostream&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Chemistry_hpp included

// ************************************************************************* //
