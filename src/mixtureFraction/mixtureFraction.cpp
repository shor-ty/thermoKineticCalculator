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

//- user def. headers
#include "MixtureFraction.hpp"

//- standard constructor
MixtureFraction::MixtureFraction() :
    MW_(0),
    Y_(0),
    Yf_(0),
    Yo_(0),
    X_(0),
    Xf_(0),
    Xo_(0),
    Con_(0),
    fuel_(false),
    oxidizer_(false)
{
}

//- one argument constructor
MixtureFraction::MixtureFraction(normalString name) : name_(name)
{
}

//- destructor
MixtureFraction::~MixtureFraction()
{
}

//- set names
void MixtureFraction::setName(const normalString& nameOfMixtureFraction)
{
    name_ = nameOfMixtureFraction;
}

//- set molecular weight
void MixtureFraction::setMW(const scalar MW)
{
    MW_ = MW;
}

//- set fuel
void MixtureFraction::setFuel()
{
    fuel_ = true;
}

//- set oxidizer
void MixtureFraction::setOxidizer()
{
    oxidizer_ = true;
}

//- return fuel
const bool MixtureFraction::fuel() const
{
    return fuel_;
}

//- return oxidizer
const bool MixtureFraction::oxidizer() const
{
    return oxidizer_;
}

//- mol fraction

    //- set mol fraction X 0 < Z 1
    void MixtureFraction::setX(const scalar& molFraction)
    {
        X_ = molFraction;
    }

    //- set mol fraction X in pure fuel
    void MixtureFraction::setXf(const scalar& molFraction)
    {
        Xf_ = molFraction;
    }

    //- set mol fraction X in pure oxidizer
    void MixtureFraction::setXo(const scalar& molFraction)
    {
        Xo_ = molFraction;
    }

    //- return mol fraction X 0 < Z < 1
    const scalar MixtureFraction::X() const
    {
        return X_;
    }

    //- return mol fraction X in pure fuel
    const scalar MixtureFraction::Xf() const
    {
        return Xf_;
    }

    //- return mol fraction X in pure oxidizer
    const scalar MixtureFraction::Xo() const
    {
        return Xo_;
    }


//- mass fraction

    //- set mass fraction Y 0 < Z < 1
    void MixtureFraction::setY(const scalar& massFraction)
    {
        Y_ = massFraction;
    }

    //- set mass fraction Y in pure fuel
    void MixtureFraction::setYf(const scalar& massFraction)
    {
        Yf_ = massFraction;
    }

    //- set mass fraction Y in pure oxidizer
    void MixtureFraction::setYo(const scalar& massFraction)
    {
        Yo_ = massFraction;
    }

    //- return mass fraction Y 0 < Z < 1
    const scalar MixtureFraction::Y() const
    {
        return Y_;
    }

    //- return mass fraction Y in pure fuel
    const scalar MixtureFraction::Yf() const
    {
        return Yf_;
    }

    //- return mass fraction Y in pure oxidizer
    const scalar MixtureFraction::Yo() const
    {
        return Yo_;
    }


//- set temperature of fuel
void MixtureFraction::setTf(const scalar& Tf)
{
    Tf_ = Tf;
}

//- set temperature of oxidizer
void MixtureFraction::setTo(const scalar& To)
{
    To_ = To;
}

//- get temperature of fuel
const scalar MixtureFraction::Tf() const
{
    return Tf_;
}

//- get temperature of oxidizer
const scalar MixtureFraction::To() const
{
    return To_;
}


const scalar MixtureFraction::MW() const
{
    return MW_;
}

//- return name
const normalString MixtureFraction::name() const
{
    return name_;
}
