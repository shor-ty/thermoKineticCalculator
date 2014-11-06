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
#include <iostream>
//- user def. headers
#include "thermodynamicProperties.hpp"

//- constructor
ThermodynamicProperties::ThermodynamicProperties()
{
}

//- destructor
ThermodynamicProperties::~ThermodynamicProperties()
{
}

//- functions
//-----------

    //- set low temperature
    void ThermodynamicProperties::setLowTemperature( const int& lowT )
    {
        T_low = lowT;
    }

    //- set mid temperature
    void ThermodynamicProperties::setMidTemperature( const int& midT )
    {
        T_mid = midT;
    }

    //- set high temperature
    void ThermodynamicProperties::setHighTemperature( const int& highT)
    {
        T_high = highT;
    }

    //- set the first line of coefficents (NASA)
    void ThermodynamicProperties::setFirstLineOfCoefficients( const double& a_1, const double& a_2, const double& a_3, const double& a_4, const double& a_5)
    {
        a1 = a_1;
        a2 = a_2;
        a3 = a_3;
        a4 = a_4;
        a5 = a_5;
    }

    //- set the second line of coefficients (NASA)
    void ThermodynamicProperties::setSecondLineOfCoefficients( const double& a_6, const double& a_7, const double& b_1, const double& b_2, const double& b_3)
    {
        a6 = a_6;
        a7 = a_7;
        b1 = b_1;
        b2 = b_2;
        b3 = b_3;
    }

    //- set the third line of coefficients (NASA)
    void ThermodynamicProperties::setThirdLineOfCoefficients( const double& b_4, const double& b_5, const double& b_6, const double& b_7)
    {
        b4 = b_4;
        b5 = b_5;
        b6 = b_6;
        b7 = b_7;
    }


//- functions
//-----------

    //- get low temperature
    int ThermodynamicProperties::getLowTemperature() const
    {
        return T_low;
    }

    //- get mid temperature
    int ThermodynamicProperties::getMidTemperature() const
    {
        return T_mid;
    }

    //- get high temperature
    int ThermodynamicProperties::getHighTemperature() const
    {
        return T_high;
    }

    //- get high temperature coefficents
    void ThermodynamicProperties::getHighTemperatureCoefficents() const
    {
        std::cout << "a1: " << a1 << "\n";
        std::cout << "a2: " << a2 << "\n";
        std::cout << "a3: " << a3 << "\n";
        std::cout << "a4: " << a4 << "\n";
        std::cout << "a5: " << a5 << "\n";
        std::cout << "a6: " << a6 << "\n";
        std::cout << "a7: " << a7 << "\n";
    }

    //- get low temperature coefficents
    void ThermodynamicProperties::getLowTemperatureCoefficents() const
    {
        std::cout << "b1: " << b1 << "\n";
        std::cout << "b2: " << b2 << "\n";
        std::cout << "b3: " << b3 << "\n";
        std::cout << "b4: " << b4 << "\n";
        std::cout << "b5: " << b5 << "\n";
        std::cout << "b6: " << b6 << "\n";
        std::cout << "b7: " << b7 << "\n";
    }

    //- calculate heat capacity cp in [kJ/kgK]
    double ThermodynamicProperties::calculateHeatCapacity( const double& T, const double& R, const double& Tlow, const double& Tmid, const double& Thigh )
    {
        //- range between low and mid temperature
        if ( T < Tmid && T >= Tlow )
        {
            return ((b1 + b2*T + b3*pow(T,2) + b4*pow(T,3) + b5*pow(T,4)) * R);
        }

        //- range between mid and high temperature
        else if ( T  > Tmid && T <= Thigh )
        {
            return ((a1 + a2*T + a3*pow(T,2) + a4*pow(T,3) + a5*pow(T,4)) * R);
        }

        //- temperature below low temperatur
        else if ( T < Tlow && T >= 273 )
        {
            std::cout << "\n   - Warning: NASA-Polynomial for heat capacity is not defined for temperature  T=" << T << " K\n   - Warning: Using low temperature coefficents for calculation...\n";
            return ((b1 + b2*T + b3*pow(T,2) + b4*pow(T,3) + b5*pow(T,4)) * R);
        }

        else
        {
            std::cout << "\n\n";
            std::cout << "   - ERROR: Temperature is out of the range of NASA-Polynomial\n";
            std::cout << "   -----------------------------------------------------------\n";
            std::cout << "     T_low = 273 K < T < " << Thigh << " K\n";
            std::cout << "     T is: " << T << " K\n\n";
            exit( EXIT_FAILURE );
        }
    }

    //- calculate heat capacity cp in [kJ/kgK]
    double ThermodynamicProperties::calculateEnthalpy( const double& T, const double& R )
    {
        return ((a1*T + a2*pow(T,2)/2 + a3*pow(T,3)/3 + a4*pow(T,4)/4 + a5*pow(T,5) + a6)*R);
    }

    //- calculate heat capacity cp in [kJ/kgK]
    double ThermodynamicProperties::calculateEntropy( const double& T, const double& R )
    {
        return ((a1*log(T) + a2*T + a3*pow(T,2)/2 + a4*pow(T,3)/3 + a5*pow(T,4)/4 + a7)*R);
    }

    void ThermodynamicProperties::showValues() const
    {
        std::cout << "\n - thermodynamic properties (NASA) \n\n";
        std::cout << "   + T(low):     " << T_low << "\n";
        std::cout << "   + T(mid):     " << T_mid << "\n";
        std::cout << "   + T(high):    " << T_high << "\n\n";
        std::cout << "   + Coeff for [T(low),T(mid)]" << "\n";
        std::cout << "   ----------------------------\n";
        std::cout << "      b1:     " << b1 << "\n";
        std::cout << "      b2:     " << b2 << "\n";
        std::cout << "      b3:     " << b3 << "\n";
        std::cout << "      b4:     " << b4 << "\n";
        std::cout << "      b5:     " << b5 << "\n";
        std::cout << "      b6:     " << b6 << "\n";
        std::cout << "      b7:     " << b7 << "\n\n";
        std::cout << "   + Coeff for [T(mid),T(high)]" << "\n";
        std::cout << "   ----------------------------\n";
        std::cout << "      a1:     " << a1 << "\n";
        std::cout << "      a2:     " << a2 << "\n";
        std::cout << "      a3:     " << a3 << "\n";
        std::cout << "      a4:     " << a4 << "\n";
        std::cout << "      a5:     " << a5 << "\n";
        std::cout << "      a6:     " << a6 << "\n";
        std::cout << "      a7:     " << a7 << "\n";
    }
