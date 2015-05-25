/*---------------------------------------------------------------------------*\
| =========                |                                                  |
| \\      /  A utomatic    | Holzmann-cfd                                     |
|  \\    /   F lamelet     | Version: 1.0                                     |
|   \\  /    C onstructor  | Web: www.Holzmann-cfd.de                         |
|    \\/                   |                                                  |
\*---------------------------------------------------------------------------*/
/*
�
�
�
�
�
�
�
�
�
�
\*---------------------------------------------------------------------------*/
//- system headers

//- user def. headers
#include "../definitions/typedef.hpp"
//#include "../thermodynamic/thermodynamic.hpp"

#ifndef FUNCTIONS_HPP_INCLUDED
#define FUNCTIONS_HPP_INCLUDED

//- for loop



//- functions

    //- open file
    const stringField openFile(const normalString&);

    //- read kinetic file
    void readChemKinThermo(const normalString&, const normalString&);

    //- calc species molecular weight [kg/mol]
    scalar calcMolecularWeight(const normalString&);


    //- create thermodynamic class objects
//    std::vector<Thermodynamic> createThermodynamicObjects(const stringField&);

#endif // FUNCTIONS_HPP_INCLUDED
