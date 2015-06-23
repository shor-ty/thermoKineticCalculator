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


class Chemistry
{
    public:

        //- constructor
        Chemistry();

        //- destructor
        ~Chemistry();


    public:

        //- function

            //- read the chemistry file
            void readChemkin
            (
                const normalString&
            );

            //- open file and return the content of the file
            stringField openFile
            (
                const normalString&
            );

            //- split string; delimiter is whitespace
            stringField splitString
            (
                const normalString
            );

            //- split reaction
            //  Returns:
            //  -  0: reactants
            //  -  1: products
            normalString splitReaction
            (
                const normalString&,
                const unsigned int
            );

            //- update stochiometric coefficient matrix nu
            void updateAllMatrix
            (
                const normalString&,
                const unsigned int,
                const scalarField&,
                const stringField&,
                const unsigned int&
            );

            //- increment the size of matrixes and vectors that are used
            void incrementMatrixesVectors();

            //- add new stochiometric factors for new elementar reaction
            //  + first int: reaction number
            //  + second int: stochiometric
            void set_nu();

            //- set kf
            void set_kf(const bool);

            //- set kb
            void set_kb(const bool);

            //- return (f)orward (b)ackward coefficient matrix
            unsigned int fb
            (
                const unsigned int,
                const unsigned int
            ) const;

            //- set formula of elementar reaction (string)
            void setElementarReaction(const normalString&);

            //- return formula of elementar reaction (string)
            normalString elementarReaction(const unsigned int&) const;

            //- insert new set of arrhenius coeffs
            void setArrheniusCoeffs
            (
                const scalar&,
                const scalar&,
                const scalar&
            );

            //- return arrhenius coeffs
            scalar arrheniusCoeffs
            (
                const unsigned int&,
                const unsigned int
            ) const;

            //- summary
            void summary() const;

    private:

        //- this field contains all elements that are used in the reaction
        stringField elements_;

        //- this field contains all species that are used in the reaction
        stringField species_;

        //- matrix that contains the reactants in elementar reaction
        stringMatrix reactants_;

        //- matrix that contains the products in elementar reaction
        stringMatrix products_;

        //- matrix of stochiometric coeffs
        matrix nu_;

        //- forward and backward reaction used or not
        scalarField kfkb_;

        //- matrix of arrhenius coeffs
        matrix arrheniusCoeffs_;

        //- TROE coefficient
        matrix TROECoeffs_;

        //- LINDEMANN formula (LOW)
        matrix LOWCoeffs_;

        //- matrix of THIRD BODY M (composition of species)
        stringMatrix Mcomp_;

        //- matrix of THIRD BODY M (values of species)
        matrix Mvalue_;

        //- vector of elementary reaction (as string)
        stringField elementarReaction_;

        //- vector of fall-off reactions
        scalarField nTROE_;

        //- vector of low-pressure reactions
        scalarField nLOW_;

        //- amount of elementar reactions
        unsigned int n_;

        //- duplicated reactions
        unsigned int nDuplicate_;

        //- TROE vector
        std::vector<int> TROE_;

        //- LOW vector
        std::vector<int> LOW_;

};

#endif // Chemistry_HPP
