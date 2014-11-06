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
#ifndef constantValues_hpp
#define constantValues_hpp
//- Debug (switch to true to use debug mode)
//------------------------------------------
bool debug{false};

//- Constants
//-----------

    //- Avogadro number [1/mol]
    double N_A = 6.0221412927e23;

    //- Bolzmannconstant [J/K]
    double k_B = 1.380648813e-23;

    //- Gas constant [J/(molK)]
    double Rm = N_A*k_B;


#endif // constantValues.hpp
