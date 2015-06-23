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

        //- this matrix contains all species of reaction r_x6
        stringMatrix speciesInReac_;

        //- matrix nu --> stochiometric coefficient matrix
        //  nu on product side === positiv
        //  nu on reactant side === negativ
        //
        //  Definition of matrix nu
        //      Range: r x m
        //      r: reactions
        //      m: species
        //  ------------------------------------------------------------------
        //     species1    species2    ...     species m
        //  ------------------------------------------------------------------
        //  1 | nu_1,1      nu_1,2      .          .
        //  2 | nu_2,1        .         .          .
        //  . |    .          .         .          .
        //  r |    .          .         .       nu_r,m
        //  ------------------------------------------------------------------
        matrix nu_;

        //- forward and backward reaction used or not
        //  Definition:
        //  -  1: both are used
        //  -  0: only forward is used
        scalarField kfkb_;

        //- arrhenius coefficient matrix nx3
        //  same as coefficient matrix for the stochiometric factors
        //
        //  -------------------------------------------------
        //  r   A0            n               Ea
        //  -------------------------------------------------
        //  1 |  x_i,j       x_i,j+1         x_i,j+2
        //  2 |  x_i+1,j        .               .
        //  3 |    .            .               .
        //  . |    .            .               .
        //  . |    .            .               .
        //  n |    .            .            x_i+n,j+2
        matrix arrheniusCoeffs_;

        //- TROE coefficient
        //  definition:
        //  +  0: T***
        //  +  1: T*
        //  +  2: T**
        matrix TROECoeffs_;

        //- LINDEMANN formula (LOW)
        //  definition:
        //  +  0: A0
        //  +  1: b
        //  +  2: Ea
        matrix LOWCoeffs_;

        //- THIRD BODY species M
        matrix M_;

        //- vector of elementary reaction (as string)
        stringField elementarReaction_;

        //- vector of fall-off reactions
        scalarField fO_;

        //- vector of low-pressure reactions
        scalarField lP_;

        //- amount of elementar reactions
        unsigned int r_;

        //- duplicated reactions
        unsigned int rDuplicate_;

        //- TROE vector
        std::vector<int> TROE_;

        //- LOW vector
        std::vector<int> LOW_;

};

#endif // Chemistry_HPP
