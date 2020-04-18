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
   AFC::ChemistryData

Description
    This class contains all chemistry data e.g. elements, reactions, species

SourceFiles
    chemistryData.cpp

\*---------------------------------------------------------------------------*/

#ifndef ChemistryData_hpp
#define ChemistryData_hpp

#include "definitions.hpp"
#include "thermo.hpp"
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

            //- Number of elementar reactions
            int nReac_{-1};

            //- List of elements
            List<word> elements_;

            //- List of species
            List<word> species_;

            //- Contains all educts species of reaction r
            List<wordList> educts_;

            //- Contains all product species of reaction r
            List<wordList> products_;

            //- Stochiometric factors for educts of reaction r
            mapList<word, int> nuEducts_;

            //- Stochiometric factors for products of reaction r
            mapList<word, int> nuProducts_;

            //- Number of duplicated reactions
            int nDuplicated_{0};

            //- Number of ignored reactions
            int nIgnored_{0};

            //- Forward reaction order for all reactions
            List<scalar> forwardReactionOrder_;

            //- Backward reaction order for all reactions
            List<scalar> backwardReactionOrder_;

            //- Global reaction order (forward + backward)
            List<scalar> globalReactionOrder_;

            //- Reaction rate kf (forward) for each reaction
            scalarField kf_;

            //- Reaction rate kb (backward) for each reaction
            scalarList kb_;

            //- Reaction rate constant Kc for each reaction
            scalarList Kc_;

            //- Elementar reaction
            List<string> elementarReaction_;

            //- Duplicated elementar reaction
            List<string> duplicatedElementarReaction_;

            //- Ignored elementar reaction
            List<string> ignoredElementarReaction_;

            //- List that contains all species in each reaction
            List<wordList> speciesInReaction_;

            //- Contains all reaction numbers where species i is included
            map<word, List<int>> reactionI_;

            //- Collision partner such as +M, +H2, +O2 for each reaction
            List<word> collisionPartner_;

            //- boolList for TBR
            //  + true if [M] is there
            //  + false if [M] is not included
            List<bool> TBR_;

            //- boolList for TBR Enhanced Factors
            //  + ture if [M] is modified
            //  + false if [M] represend all species
            List<bool> ENHANCE_;

            //- boolList for TBR LOW
            //  + true if low pressure modifications using Lindemann formula
            //  + false if common reaction
            List<bool> LOW_;

            //- boolList for TBR TROE (also need LOW coeffs)
            //  + true if low pressure modifications using TROE formula
            //  + false if common reaction
            List<bool> TROE_;

            //- boolList for TBR SRI
            List<bool> SRI_;

            //- boolList for forward reaction only =>
            //  This list even contain the information of irreversible reac.
            List<bool> forwardReaction_;

            //- boolList for backward reactions only <=
            //  This list even contain the information of irreversible reac.
            List<bool> backwardReaction_;

            //- Matrix of Arrhenius coeffs
            List<scalarField> arrheniusCoeffs_;

            //- Matrix of TROE coeffs
            List<scalarField> TROECoeffs_;

            //- Matrix of Arrhenius coeffs for high pressure
            List<scalarField> LOWCoeffs_;

            //- Matrix of SRI coeffs
            List<scalarField> SRICoeffs_;

            //- Enhanced coeffs for adjustment
            mapList<word, scalar> ENHANCEDCoeffs_;

            //- Omega
            scalarField omega_;

            scalar dH_;
            scalar dS_;
            scalar dG_;


        //- Thermodynamic available in chemistry file
        //bool thermo_;


    public:

        //- Constructor
        ChemistryData(const string);

        //- Destructor
        ~ChemistryData();

        // Member functions

            //- Set the thermo bool (if true the NASA polynomials are
            //  inside the chemistry file)
            //void setThermo();

            //- Return the thermo switch
            //bool thermo() const;


        // Insert functions, from ChemistryReader:: delegated

            //- Insert elements
            void elements(const word);

            //- Insert species
            void species(const word);

            //- Insert species as educt
            void educt(const word);

            //- Insert species as product
            void product(const word);

            //- Insert stochiometric coeff for educts of actual reaction
            void nuEducts(const word, const int);

            //- Insert stochiometric coeff for products of actual reaction
            void nuProducts(const word, const int);

            //- Insert elementar reaction
            void elementarReaction(const string);

            //- Insert duplicated elementar reaction
            void duplicatedElementarReaction(const string);

            //- Insert ignored elementar reaction
            void ignoredElementarReaction(const string);

            //- Set arrhenius coeffs
            void arrheniusCoeffs(const scalar, const scalar, const scalar);

            //- Set collision partner
            void collisionPartner(const word);

            //- Insert LOW coeffs
            void LOWCoeffs(const scalar, const unsigned int);

            //- Insert TROE coeffs
            void TROECoeffs(const scalar, const unsigned int);

            //- Insert SRI coeffs
            void SRICoeffs(const scalar, const unsigned int);

            //- Insert ENHANCE factors (species + value)
            void ENHANCEDCoeffs(const word, const scalar);

            //- Increment nDuplicated_
            void incrementDuplicated();

            //- Increment nIgnored
            void incrementIgnored();

            //- Increment nReac_
            void incrementReac();

            //- Build List<wordList> of species included in reaction r
            void speciesInReaction(const word);

            //- Increment vectors and matrix size
            void incrementMatrixesVectors();


        // Setter functions, from ChemistryReader:: delegated

            //- Set [M]
            void setM
            (
                const scalar&
            );

            //- Set the forward reaction boolean
            void FR(const bool);

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
            void setReacNumbers(const word, const int);

            //- Set forward reaction order
            void forwardReactionOrder();

            //- Set backward reaction order
            void backwardReactionOrder();


        // Update functions

            //- Update reaction rates kf
            void updateKf(const int, const scalar);

            //- Update reaction rates kb
            void updateKb(const int, const scalar);

            //- Update reaction rates constant Kc
            void updateKc(const int, const scalar);

            //- Update source term omega for species s
            void calculateOmega(const int, const scalar);

            //- Update source term omega for all species
            void calculateOmega(const scalarField&);

            //- Update the global reaction order list
            void updateGlobalReactionOrder();

        // Return functions


            //- Return bool

                //- Forward reaction
                bool FR(const int) const;

                //- Backward reaction
                bool BR(const int) const;

                //- Return true if elementar reaction is a TBR reaction
                bool TBR(const int);

                //- Return true if elementar reaction is a TBR reaction
                bool TBR(const int) const;

                //- Return true if elementar reaction is a LOW reaction
                bool LOW(const int) const;

                //- Return true if elementar reaction is a TROE reaction
                bool TROE(const int) const;

                //- Return true if elementar reaction is a SRI reaction
                bool SRI(const int) const;

                //- Return true if elementar reaction has enhanced factors
                bool ENHANCED(const int) const;


            //- Return M_ [g/cm^3]
            scalar M() const;

            //- Return dH
            scalar dS() const;

            //- Return dG
            scalar dG() const;

            //- Return all elements
            wordList elements() const;

            //- Return all species
            wordList species() const;

            //- Return the educt species of reaction r
            wordList educts(const int) const;

            //- Return the product species of reaction r
            wordList products(const int) const;

            //- Return amount of duplicated reactions
            unsigned int nDuplicated() const;

            //- Return amount of ignored reactions
            unsigned int nIgnored() const;

            //- Return no. of reaction
            int nReac() const;

            //- Return all elementar reactions in chemistry
            stringList elementarReaction() const;

            //- Return elementar reaction (as string)
            string elementarReaction(const int) const;

            //- Return ignored elementar reaction
            wordList ignoredElementarReaction() const;

            //- Return ignored elementar reaction (as string)
            word ignoredElementarReaction(const int) const;

            //- Return List of reaction no. of species
            List<int> reacNumbers(const word) const;

            //- Return species list for reaction r
            List<wordList> speciesInReaction() const;

            //- Return species list for reaction r
            wordList speciesInReaction(const int) const;

            //- Return stochiometric factors of educts of reaction r
            map<word, int> nuEducts(const int) const;

            //- Return stochiometric factors of educts of reaction r
            map<word, int> nuProducts(const int) const;

            //- Return the exponent factor for Keq calculation
            scalar exponent(const int) const;

            //- Return if forward reaction possible
            bool forwardReaction(const int) const;

            //- Return if backward reaction possible
            bool backwardReaction(const int) const;

            //- Return arrhenius coeffs for reaction no.
            //  [0] -> pre-exponent [units depend on equation]
            //  [1] -> temperature exponent [-]
            //  [2] -> activation energy [cal/mol]
            scalarList arrheniusCoeffs(const int) const;

            //- Return the collision number of reaction no.
            word collisionPartner(const int) const;

            //- Return arrhenius coeffs for high pressure for reaction no.
            scalarList LOWCoeffs(const int) const;

            //- Return TROE coeffs
            scalarList TROECoeffs(const int) const;

            //- Return SRI coeffs
            scalarList SRICoeffs(const int) const;

            //- Return ENHANCED factors (species + value) of reac no.
            map<word, scalar> ENHANCEDCoeffs (const int) const;

            //- Return reaction rate kf for reaction no.
            scalar kf(const int) const;

            //- Return reaction rates kf
            scalarList kf() const;

            //- Return reaction rate kb for reaction no.
            scalar kb(const int) const;

            //- Return reaction rates kb
            scalarList kb() const;

            //- Return reaction rate constant Kc for reaction no.
            scalar Kc(const int) const;

            //- Return reaction rate constant Kc
            scalarList Kc() const;

            //- Return omega of species s
            scalar omega(const int) const;

            //- Return omega field
            scalarField omega() const;

            //- Return forward reaction order of reaction r
            scalar forwardReactionOrder(const int) const;

            //- Return backward reaction order of reaction r
            scalar backwardReactionOrder(const int) const;

            //- Return global reaction order of reaction r
            scalar globalReactionOrder(const int) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ChemistryData_hpp included

// ************************************************************************* //
