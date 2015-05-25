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
#include "thermodynamic.hpp"

//- constructor
Thermodynamic::Thermodynamic()
{
}

//- destructor
Thermodynamic::~Thermodynamic()
{
}

//- functions
//-----------


    //- set phase
    void Thermodynamic::setPhase( const normalString& phase_ )
    {
        phase = phase_;
    }

    //- set low temperature
    void Thermodynamic::setPolyTemperature
    (
        const scalar& lowTemp_,
        const scalar& comTemp_,
        const scalar& higTemp_
    )
    {
        lowTemp = lowTemp_;
        comTemp = comTemp_;
        higTemp = higTemp_;
    }


    //- set the first line of coefficents (NASA)
//    void Thermodynamic::setNASACoefficients( const double& a_1, const double& a_2, const double& a_3, const double& a_4, const double& a_5)
//    {
//
//    }

//- functions
//-----------

    //- get low temperature
//    int Thermodynamic::lowTemperature() const
//    {
//        return T_low;
//    }
//
//    //- get mid temperature
//    int Thermodynamic::midTemperature() const
//    {
//        return T_mid;
//    }
//
//    //- get high temperature
//    int Thermodynamic::highTemperature() const
//    {
//        return T_high;
//    }
//
    //- get high temperature coefficents
    void Thermodynamic::setPolyCoeffs(const scalarField& polyCoeffs_)
    {
        polyCoeffs = polyCoeffs_;
    }
//
//    //- get low temperature coefficents
//    void Thermodynamic::lowTemperatureCoefficents() const
//    {
//        std::cout << "b1: " << b1 << "\n";
//        std::cout << "b2: " << b2 << "\n";
//        std::cout << "b3: " << b3 << "\n";
//        std::cout << "b4: " << b4 << "\n";
//        std::cout << "b5: " << b5 << "\n";
//        std::cout << "b6: " << b6 << "\n";
//        std::cout << "b7: " << b7 << "\n";
//    }
//
    //- calculate heat capacity cp in [kJ/kgK]
    scalar Thermodynamic::calculateHeatCapacity(const scalar& T_)
    {
        int i{0};

        //- get information of temperature range
        whichPolyCoeffs(T_, i);

        // return heat capacity [J/(molK)]
        return ((   polyCoeffs[i]
                  + polyCoeffs[i+1]*T_
                  + polyCoeffs[i+2]*pow(T_,2)
                  + polyCoeffs[i+3]*pow(T_,3)
                  + polyCoeffs[i+4]*pow(T_,4) )
                  * R);
    }

    //- calculate enthalpy
    scalar Thermodynamic::calculateEnthalpy(const scalar& T_)
    {
        int i{0};

        //- get information of temperature range
        whichPolyCoeffs(T_, i);

        // return enthalpy [J/(mol)]
        return ((   polyCoeffs[i]
                  + polyCoeffs[i+1]*T_/2
                  + polyCoeffs[i+2]*pow(T_,2)/3
                  + polyCoeffs[i+3]*pow(T_,3)/4
                  + polyCoeffs[i+4]*pow(T_,4)
                  + polyCoeffs[i+5]/T_)
                  *T_*R);
    }

    //- calculate entropy
    scalar Thermodynamic::calculateEntropy(const scalar& T_)
    {
        int i{0};

        //- get information of temperature range
        whichPolyCoeffs(T_, i);

        // return entropy [J/(molK)]
        return ((   polyCoeffs[i]*log(T_)
                  + polyCoeffs[i+1]*T_
                  + polyCoeffs[i+2]*pow(T_,2)/2
                  + polyCoeffs[i+3]*pow(T_,3)/3
                  + polyCoeffs[i+4]*pow(T_,4)/4
                  + polyCoeffs[i+6])
                  *R);
    }


    //- set thermodynamic bool to true
    void Thermodynamic::thermodynamicTrue()
    {
        thermoBool = true;
    }


    //- return thermodynamic status
    bool Thermodynamic::thermodynamicStatus() const
    {
        return thermoBool;
    }

    //- check temperature range
    void Thermodynamic::whichPolyCoeffs(const scalar& T_, int& i)
    {
        //- FIXME EXTEND TO WHICH SPECIES
        if(T_ < lowTemp)
        {
            std::cerr << "    - Warning: NASA-Polynomial is not defined for temperature  T=" << T_ << " K\n    \
                              - Warning: Using low temperature coefficents for calculation...\n";
            i = 7;
        }
        else if(T_ >= lowTemp && T_  <= comTemp)
        {
            i = 7;
        }
        else if(T_ > comTemp && T_  <= higTemp)
        {
            i = 0;
        }
        else
        {
            std::cerr << "    - Warning: NASA-Polynomial is not defined for temperature  T=" << T_ << " K\n    \
                              - Warning: Using high temperature coefficents for calculation...\n";
            i = 0;
        }
    }
//    void Thermodynamic::showValues() const
//    {
//        std::cout << "\n - thermodynamic  (NASA) \n\n";
//        std::cout << "   + T(low):     " << T_low << "\n";
//        std::cout << "   + T(mid):     " << T_mid << "\n";
//        std::cout << "   + T(high):    " << T_high << "\n\n";
//        std::cout << "   + Coeff for [T(low),T(mid)]" << "\n";
//        std::cout << "   ----------------------------\n";
//        std::cout << "      b1:     " << b1 << "\n";
//        std::cout << "      b2:     " << b2 << "\n";
//        std::cout << "      b3:     " << b3 << "\n";
//        std::cout << "      b4:     " << b4 << "\n";
//        std::cout << "      b5:     " << b5 << "\n";
//        std::cout << "      b6:     " << b6 << "\n";
//        std::cout << "      b7:     " << b7 << "\n\n";
//        std::cout << "   + Coeff for [T(mid),T(high)]" << "\n";
//        std::cout << "   ----------------------------\n";
//        std::cout << "      a1:     " << a1 << "\n";
//        std::cout << "      a2:     " << a2 << "\n";
//        std::cout << "      a3:     " << a3 << "\n";
//        std::cout << "      a4:     " << a4 << "\n";
//        std::cout << "      a5:     " << a5 << "\n";
//        std::cout << "      a6:     " << a6 << "\n";
//        std::cout << "      a7:     " << a7 << "\n";
//    }
