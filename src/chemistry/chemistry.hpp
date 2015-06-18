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

            //- add new stochiometric factors for new elementar reaction
            //  + first int: reaction number
            //  + second int: stochiometric
            void addNu(const int&, const int&);

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

            //- increment the amount of reactions
            void increment_r();

            //- decrement the amount of reactions
            //  used for LOW TROE
            void decrement_r();

            //- return the amount of elementar reactions
            int r() const;

            //- set TROE vector
            void setTROE(const bool);

            //- set LOW vector
            void setLOW(const bool);

            //- return TROE_ coefficient vector
            unsigned int TROE(const unsigned int) const;

            //- return LOW_ coefficient vector
            unsigned int LOW(const unsigned int) const;

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



    private:

        //- stochiometric factors nu
        //  as a coefficient matrix
        //  maximal third body reaction means 3 reactants 3 products = 6
        //  matrix is nx6 := n is the amount of elementar Chemistry
        //
        //  ------------------------------------------------------------------
        //  r   elem1    ...    elem3          elem4        ... elem6
        //  ------------------------------------------------------------------
        //  1 | nu_i,j   ...    nu_i,j+2        nu_i,j+3    ... nu_i,j+5
        //  2 | nu_i+1,j ...    nu_i+1,j+2      nu_i+1,j+3  ... nu_i+1,j+5
        //  3 | nu_i+2,j ...        .               .               .
        //  . |    .                .               .               .
        //  . |    .                .               .               .
        //  n |    .                .               .           nu_i+n,j+5
        //  ------------------------------------------------------------------
        matrix nu_;

        //- forward and backward reaction used or not
        //  matrix
        //
        //  r     kf     kb
        //  -----------------
        //  1  |  1       0
        //  2  |  1       1
        //  .  |  .       .
        //  n  |  .       .
        //  -----------------
        matrix fb_;

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

        //- matrix of elementary reaction (as string)
        //  only for storing - not needed afterwards
        public:
        std::vector<normalString> elementarReaction_;

        //- TROE vector
        std::vector<int> TROE_;

        //- LOW vector
        std::vector<int> LOW_;

        //- amount of elementar Chemistry
        unsigned int r_;
};

#endif // Chemistry_HPP
