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
#ifndef ThermodynamicProperties_hpp
#define ThermodynamicProperties_hpp
//- system headers
#include <math.h>
//- user def. headers

class ThermodynamicProperties
{
    public:

        //- constructor
        ThermodynamicProperties();

        //- destructor
        ~ThermodynamicProperties();


    public:

    //- functions
    //-----------

        //- set low temperatur
        void setLowTemperature( const int & );

        //- get low temperatur
        int getLowTemperature() const;

        //- set mid temperature
        void setMidTemperature( const int& );

        //- get mid temperature
        int getMidTemperature() const;

        //- set high temperature
        void setHighTemperature( const int& );

        //- get high temperatur
        int getHighTemperature() const;

        //- show values
        void showValues() const;

        //- set first line of coefficents
        void setFirstLineOfCoefficients( const double&, const double&, const double&, const double&, const double& );

        //- set second line of coefficents
        void setSecondLineOfCoefficients( const double&, const double&, const double&, const double&, const double& );

        //- set third line of coefficents
        void setThirdLineOfCoefficients( const double&, const double&, const double&, const double& );

        //- get all low temperature coefficents (NASA-Polynomials)
        void getLowTemperatureCoefficents() const;

        //- get all high temperature coefficents (NASA-Polynomials)
        void getHighTemperatureCoefficents() const;

        //- calculate heat capacity cp [kJ/kgK]
        double calculateHeatCapacity( const double&, const double&, const double&, const double&, const double& );

        //- calculate heat capacity H []
        double calculateEnthalpy( const double&, const double& );

        //- calculate heat capacity S []
        double calculateEntropy( const double&, const double& );


    private:

    //- variables
    //-----------

        //- low temperature
        int T_low{0};

        //- mid temperature
        int T_mid{0};

        //- high temperature
        int T_high{0};

        //- high temperature coefficients (NASA-P)
        double a1{0.}, a2{0.}, a3{0.}, a4{0.}, a5{0.}, a6{0.}, a7{0.};

        //- low temperature coefficient (NASA-P)
        double b1{0.}, b2{0.}, b3{0.}, b4{0.}, b5{0.}, b6{0.}, b7{0.};
};


#endif // ThermodynamicProperties_hpp
