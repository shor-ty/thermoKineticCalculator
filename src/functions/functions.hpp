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
//- system headers

//- user def. headers
#include "../definitions/typedef.hpp"
#include "../mixtureFraction/mixtureFraction.hpp"


#ifndef FUNCTIONS_HPP_INCLUDED
#define FUNCTIONS_HPP_INCLUDED

//- for loop



//- functions

    //- read mixtureFractionPoints
    const unsigned int mixtureFractionPoints
    (
        const normalString&
    );

    //- read scalar dissipation rates
    const scalarField scalarDissipationRates
    (
        const normalString&
    );

    //- read mixtureFractionPoints




    //- calculate stochiometric mixture fraction
//    scalar stochiometricMF
//    (
//        const std::vector<Species>&
//    );
//
//    //- calculate adiabatic enthalpy of fuel and oxidizer
//    //  + 0 mean fuel
//    //  + 1 mean oxidizer
//    //  [J/kg]
//    scalar adiabaticEnthalpy
//    (
//        const std::vector<Species>&,
//        const int
//    );
//
//    //- calculate adiabatic flame temperature at stochiometric
//    //  mixture in Kelvin
//    scalar adiabateFlameTemperature
//    (
//        const scalar&,
//        const std::vector<Species>&
//    );
//
//    //- discret points of mixture fraction points Z [-]
//    scalarField discretZ
//    (
//        const normalString&,
//        const scalar&
//    );
//
//    //- discrete points of scalar dissipation rate [Hz]
//    scalarField discretChi
//    (
//        const normalString&
//    );
//
//    //- calculate the linear distribution of the mass fraction of
//    //  fuel and oxidizer at the discrete points
//    void initializeY
//    (
//        const scalarField&,
//        std::vector<std::vector<Species> >&
//    );
//
//    //- remove white space
//    void removeSpace
//    (
//        normalString&
//    );
//
//    //- calc species molecular weight [kg/mol]
//    scalar calcMolecularWeight
//    (
//        const normalString&
//    );
//
//    //- summary of data
//    void summary
//    (
//        const scalar&,
//        const scalar&,
//        const scalar&,
//        const scalarField&
//        //const scalar&,
//       // const std::vector<Species>&
//    );

//---------------------------------------------------------------------------*/

#endif // FUNCTIONS_HPP_INCLUDED

//---------------------------------------------------------------------------*/


