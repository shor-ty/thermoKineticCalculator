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
:
    thermoBool(false)
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


//- functions
//-----------

    //- set NASA coefficents
    void Thermodynamic::setPolyCoeffs(const scalarField& polyCoeffs_)
    {
        polyCoeffs = polyCoeffs_;
    }

    //- calculate heat capacity cp [J/(molK)]
    scalar Thermodynamic::cp(const scalar T_) const
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

    //- calculate enthalpy [J/mol]
    scalar Thermodynamic::h(const scalar T_) const
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

    //- calculate entropy [J/(molK)]
    scalar Thermodynamic::s(const scalar T_) const
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
    void Thermodynamic::whichPolyCoeffs(const scalar& T_, int& i) const
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
