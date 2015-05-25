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
//- system headers
#include <iostream>
//- user def. headers
#include "elements.hpp"

Elements::Elements() :
elementName(6), elementAM(6)
{
    //- Initialize element names

        elementName[0] = "H";
        elementName[1] = "O";
        elementName[2] = "C";
        elementName[3] = "N";
        elementName[4] = "AR";
        elementName[5] = "HE";

    //- Initialize atomic mass

        elementAM[0] = 1.008;       // H
        elementAM[1] = 15.999;      // O
        elementAM[2] = 12.011;      // C
        elementAM[3] = 14.007;      // N
        elementAM[4] = 39.948;      // AR
        elementAM[5] = 4.002602;    // HE
}

Elements::~Elements()
{
}

scalar Elements::atomicWeight(const normalString element)
{
    //- ID
    int id{-1};

    //- Search element
    for (unsigned int i=0; i<elementName.size(); i++)
    {
        if (element == elementName[i])
        {
            id = i;
        }
    }

    //- if id still -1 abort
    if (id == -1)
    {
        std::cerr << "\n ++ ERROR: Element " << element << " is not in the database; atomic weight is missing...";
        std::cerr << "\n ++ Please email me: Tobias.Holzmann@Holzmann-cfd.de";
        std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__ << std::endl;
        std::terminate();
    }

    return elementAM[id];
}
