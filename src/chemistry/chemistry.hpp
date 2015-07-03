/*---------------------------------------------------------------------------*\
| =========                |                                                  |
| \\      /  A utomatic    | Holzmann-cfd                                     |
|  \\    /   F lamelet     | Version: 1.0                                     |
|   \\  /    C onstructor  | Web: www.Holzmann-cfd.de                         |
|    \\/                   |                                                  |
\*---------------------------------------------------------------------------*/
/*
»
»
»
»
»
»
»
»
»
»
\*---------------------------------------------------------------------------*/
#ifndef Chemistry_hpp
#define Chemistry_hpp
//- system headers


//- user def. headers
#include "../definitions/typedef.hpp"
#include "../thermodynamic/thermodynamic.hpp"
#include "../transport/transport.hpp"


class Chemistry
:
    public Thermodynamic
{
    public:

        //- constructor
        Chemistry();

        //- destructor
        ~Chemistry();


    public:

        //- functions

            //- read the chemistry file
            void readChemkin
            (
                const normalString&
            );

            //- increment the size of matrixes and vectors that are used
            void incrementMatrixesVectors();

            //- save elementar reaction as re-arranged string
            void elementarReaction
            (
                const normalString&
            );

            //- set arrhenius coeffs
            void arrheniusCoeffs
            (
                const normalString&
            );

            //- update ENHANCED factors matrix
            void enhancedFactors
            (
                const normalString&
            );

            //- update LOW pressure arrhenius coeffs matrix
            void LOWCoeffs
            (
                const normalString&
            );

            //- update TROE coeffs matrix
            void TROECoeffs
            (
                const normalString&
            );

            //- update SRI coeffs matrix
            void SRICoeffs
            (
                const normalString&
            );

            //- update backward reaction vector
            void backwardReaction();

            //- update reactants and products matrix (species)
            void reactantsAndProducts();

            //- update stochiometric matrix (pre-processing)
            void updateStochiometricMatrix
            (
                const stringField&,
                const unsigned int,
                const unsigned int&
            );

            //- update stochiometric matrix (store data)
            void updateStochiometricMatrix
            (
                const normalString&,
                const scalar,
                const unsigned int,
                const unsigned int&
            );

            //- remove THIRD BODY from reaction
            void removeThirdBody
            (
                normalString&,
                const unsigned int&
            );


            //- summary
            void summary() const;

            //- read thermo data from chemistry file
            void readChemKinThermo
            (
                const normalString&
            );

            //- create reaction rate matrix
            void createReactionRateMatrix();

            //- check thermo
            void checkThermo
            (
                const normalString&
            ) const;

            //- check thermo
            void checkTrans
            (
                const normalString&
            ) const;

            //- return thermodynamic_ value
            bool thermo() const;

            //- return species field
            stringField species() const;



    private:

        //- vector of all elements that are used in the chemistry
        stringField elements_;

        //- vector of all species that are used in the chemistry
        stringField species_;

        //- matrix of stochiometric coeffs
        matrix nu_;

        //- matrix of arrhenius coeffs
        matrix arrheniusCoeffs_;

        //- vector of elementary reaction (as string)
        stringField elementarReaction_;

        //- amount of elementar reactions
        int n_;

        //- duplicated reactions
        unsigned int nDuplicate_;

        //- vector of backward reaction
        std::vector<bool> kb_;


        //- vector for THIRD BODY REACTION
        std::vector<bool> TBR_;


        //- SRI THIRD BODY

            //- vector for SRI
            std::vector<bool> SRI_;

            //- SRI coeffs
            matrix SRICoeffs_;


        //- TROE THIRD BODY

            //- vector for TROE
            std::vector<bool> TROE_;

            //- TROE coefficient
            matrix TROECoeffs_;


        //- LOW THIRD BODY

            //- vector for low pressure
            std::vector<bool> LOW_;

            //- LOW pressure arrhenius coeffs
            matrix LOWCoeffs_;


        //- ENHANCED THIRD BODY

            //- vector for THIRD BODY REACTION for ENHANCE
            std::vector<bool> ENHANCE_;

            //- matrix of THIRD BODY M (composition of species)
            stringMatrix Mcomp_;

            //- matrix of THIRD BODY M (values of species)
            matrix Mvalue_;


        //- reaction rate matrix (this matrix contain the
        //  reaction no. that influence the species
        matrix reactionRateMatrix_;

        //- bool if thermodynamic is included in the chemkin file
        bool themodynamic_;

};
// ************************************************************************* //
#endif // Chemistry_HPP
