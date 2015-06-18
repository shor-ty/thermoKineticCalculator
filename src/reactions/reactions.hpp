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
#ifndef Reactions_hpp
#define Reactions_hpp
//- system headers


//- user def. headers
#include "../definitions/typedef.hpp"


class Reactions
{
    public:

        //- constructor
        Reactions();

        //- destructor
        ~Reactions();

    public:

        //- function

            //- set kf
            void set_kf();

            //- set kb
            void set_kb();

            //- return status of kf
            bool kf() const;

            //- return status of kb
            bool kb() const;

            //- set elementar reaction
            void setElementarReaction(const normalString&);

            //- return elementar reaction
            normalString elementarReaction() const;

            //- add new stochiometric factors for new elementar reaction
            //  + first int: reaction number
            //  + second int: stochiometric
            void addNu(const int&, const int&);

            //- create matrixes
            void createMatrix(const normalString);


            //- set arrhenius coefficients
            void setArrheniusCoeffs
            (
                const scalar&,
                const scalar&,
                const scalar&
            );

            //- set TROE coefficients
            void setTROECoeffs
            (
                const scalar&,
                const scalar&,
                const scalar&,
                const scalar&
            );

            //- set LOW
            void setLOW();

            //- return status of LOW
            bool statusLOW() const;

            //- set TROE
            void setTROE();

            //- return status of TROE
            bool statusTROE() const;

            //- return TROE coefficient a
            scalar TROE_a() const;

            //- return TROE coefficient T*
            scalar TROE_Ts() const;

            //- return TROE coefficient T**
            scalar TROE_Tss() const;

            //- return TROE coefficient T***
            scalar TROE_Tsss() const;

            //- return activation energy
            scalar arrhenius_Ea() const;

            //- return prä-exponential factor (frequency factor)
            scalar arrhenius_A() const;

            //- return temperature correctur factor
            scalar arrhenius_b() const;

            //- return LOW activation energy
            scalar arrhenius_Ea_LOW() const;

            //- return LOW prä-exponential factor (frequency factor)
            scalar arrhenius_A_LOW() const;

            //- return LOW temperature correctur factor
            scalar arrhenius_b_LOW() const;


    private:

        //- stochiometric factors nu
        //  as a matrix
        //  maximal third body reaction means 3 reactants 3 products = 6
        //  matrix is nx6 := n is the amount of elementar reactions
        //
        //     ---------------------------------------------------------------
        //      elem1    ...    elem3          elem4        ... elem6
        //     ---------------------------------------------------------------
        //  1 | nu_i,j   ...    nu_i,j+2        nu_i,j+3    ... nu_i,j+5
        //  2 | nu_i+1,j ...    nu_i+1,j+2      nu_i+1,j+3  ... nu_i+1,j+5
        //  3 | nu_i+2,j ...        .               .               .
        //  . |    .                .               .               .
        //  . |    .                .               .               .
        //  n |    .                .               .           nu_i+n,j+5
        //     ---------------------------------------------------------------
        matrix nu_;

        //- amount of elementar reactions
        int r_;

        //- elementar reaction
        normalString elementarReaction_;

        //- bool forward reaction
        bool kf_;

        //- bool backward reaction
        bool kb_;

        //- bool LOW
        bool LOW;

            //- LOW coefficient

                //- activation energy
                scalar Ea_LOW;

                //- prä-exponential factor (frequency factor)
                scalar A_LOW;

                //- temperature correctur factor
                scalar b_LOW;


        //- bool TROE
        bool TROE;

            //- TROE coefficient

                //- correction coefficient a
                scalar a;

                //- correction coefficient T*
                scalar Ts;

                //- correction coefficient T**
                scalar Tss;

                //- correction coefficient T***
                scalar Tsss;


        //- arrhenius constants
        //  + if LOW available these represent for high-pressure

            //- activation energy
            scalar Ea;

            //- prä-exponential factor (frequency factor)
            scalar A;

            //- temperature correctur factor
            scalar b;

        //- LOW extension

        //- TROE extention
};

#endif // REACTIONS_HPP
