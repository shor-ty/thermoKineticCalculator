/*---------------------------------------------------------------------------*\
| =========                |                                                  |
| \\      /  A utomatic    | Holzmann-cfd                                     |
|  \\    /   F lamelet     | Version: 1.0                                     |
|   \\  /    C onstructor  | Web: www.Holzmann-cfd.de                         |
|    \\/                   |                                                  |
\*---------------------------------------------------------------------------*/
/*
»
» Description:
»   This class contains information about the elements:
»       + Name
»       + Molecular weight
»
»
» Used:
»   For calculating the net rate of Elements
»
\*---------------------------------------------------------------------------*/
#ifndef Elements_HPP
#define Elements_HPP
//- system headers
#include <vector>
#include <string>

//- user def. headers
#include "../definitions/typedef.hpp"


class Elements
{
    public:

        //- constructor
        Elements();

        //- destructor
        ~Elements();


    public:

        //- functions

            //- return atomic mass of element [kg/mol]
            scalar atomicWeight
            (
                const normalString&,
                const scalar&
            );


    private:

        //- database array:

            //- name of element
            stringField elementName;

            //- atomic mass of element
            scalarField elementAM;
};

#endif // Elements_HPP
