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
#ifndef Thermodynamic_hpp
#define Thermodynamic_hpp
//- system headers
#include <math.h>
//- user def. headers
#include "../definitions/typedef.hpp"


class Thermodynamic
{
    public:

        //- constructor
        Thermodynamic();

        //- destructor
        ~Thermodynamic();


    public:

    //- functions
    //-----------


        //- set phase
        void setPhase( const normalString& );

        //- set temperature bounds
        void setPolyTemperature( const scalar&, const scalar&, const scalar& );


//
//        //- get low temperatur
//        int lowTemperature() const;
//
//        //- set mid temperature
//        void setMidTemperature( const int& );
//
//        //- get mid temperature
//        int midTemperature() const;
//
//        //- set high temperature
//        void setHighTemperature( const int& );
//
//        //- get high temperatur
//        int highTemperature() const;
//
//        //- show values
//        void showValues() const;
//
//        //- set polyCoeffs
          void setPolyCoeffs( const scalarField& );
//
//        //- set second line of coefficents
//        void setSecondLineOfCoefficients( const double&, const double&, const double&, const double&, const double& );
//
//        //- set third line of coefficents
//        void setThirdLineOfCoefficients( const double&, const double&, const double&, const double& );
//
//        //- get all low temperature coefficents (NASA-Polynomials)
//        void lowTemperatureCoefficents() const;
//
//        //- get all high temperature coefficents (NASA-Polynomials)
//        void highTemperatureCoefficents() const;
//
        //- return heat capacity cp [J/(molK)]
        scalar cp( const scalar ) const;

        //- return enthalpy [J/(mol)]
        scalar h( const scalar ) const;

        //- return ref enthalpy [J/mol]
        scalar h0() const;

        //- calculate entropie S [J/(molK)]
        scalar s( const scalar ) const;

        //- calculate free gibbs energy
        scalar G( const scalar ) const ;

        //- set thermodynamic status
        void thermodynamicTrue();

        //- return thermodynamic status
        bool thermodynamicStatus() const;

        //- define the
        void whichPolyCoeffs( const scalar&, int& ) const;


    private:

    //- variables
    //-----------

        //- phase
        normalString phase;

        //- themodynamic values set
        bool thermoBool;

        //- temperature bounds
        scalar lowTemp;
        scalar comTemp;
        scalar higTemp;

        //- polyCoeffs
        scalarField polyCoeffs;
//
//        //- mid temperature
//        int T_mid{0};
//
//        //- high temperature
//        int T_high{0};
//
//        //- high temperature coefficients (NASA-P)
//        double a1{0.}, a2{0.}, a3{0.}, a4{0.}, a5{0.}, a6{0.}, a7{0.};
//
//        //- low temperature coefficient (NASA-P)
//        double b1{0.}, b2{0.}, b3{0.}, b4{0.}, b5{0.}, b6{0.}, b7{0.};
};


#endif // Thermodynamic_hpp
