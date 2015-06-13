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
#include "species.hpp"

//- standard constructor
Species::Species() :
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
Species::Species(normalString name) : name_(name)
{
}

//- destructor
Species::~Species()
{
}

//- set names
void Species::setName(const normalString& nameOfSpecies)
{
    name_ = nameOfSpecies;
}

//- set molecular weight
void Species::setMW(const scalar MW)
{
    MW_ = MW;
}

//- set fuel
void Species::setFuel()
{
    fuel_ = true;
}

//- set oxidizer
void Species::setOxidizer()
{
    oxidizer_ = true;
}

//- return fuel
const bool Species::fuel() const
{
    return fuel_;
}

//- return oxidizer
const bool Species::oxidizer() const
{
    return oxidizer_;
}

//- mol fraction

    //- set mol fraction X 0 < Z 1
    void Species::setX(const scalar& molFraction)
    {
        X_ = molFraction;
    }

    //- set mol fraction X in pure fuel
    void Species::setXf(const scalar& molFraction)
    {
        Xf_ = molFraction;
    }

    //- set mol fraction X in pure oxidizer
    void Species::setXo(const scalar& molFraction)
    {
        Xo_ = molFraction;
    }

    //- return mol fraction X 0 < Z < 1
    const scalar Species::X() const
    {
        return X_;
    }

    //- return mol fraction X in pure fuel
    const scalar Species::Xf() const
    {
        return Xf_;
    }

    //- return mol fraction X in pure oxidizer
    const scalar Species::Xo() const
    {
        return Xo_;
    }


//- mass fraction

    //- set mass fraction Y 0 < Z < 1
    void Species::setY(const scalar& massFraction)
    {
        Y_ = massFraction;
    }

    //- set mass fraction Y in pure fuel
    void Species::setYf(const scalar& massFraction)
    {
        Yf_ = massFraction;
    }

    //- set mass fraction Y in pure oxidizer
    void Species::setYo(const scalar& massFraction)
    {
        Yo_ = massFraction;
    }

    //- return mass fraction Y 0 < Z < 1
    const scalar Species::Y() const
    {
        return Y_;
    }

    //- return mass fraction Y in pure fuel
    const scalar Species::Yf() const
    {
        return Yf_;
    }

    //- return mass fraction Y in pure oxidizer
    const scalar Species::Yo() const
    {
        return Yo_;
    }


//- set temperature of fuel
void Species::setTf(const scalar& Tf)
{
    Tf_ = Tf;
}

//- set temperature of oxidizer
void Species::setTo(const scalar& To)
{
    To_ = To;
}

//- get temperature of fuel
const scalar Species::Tf() const
{
    return Tf_;
}

//- get temperature of oxidizer
const scalar Species::To() const
{
    return To_;
}


const scalar Species::MW() const
{
    return MW_;
}

//- return name
const normalString Species::name() const
{
    return name_;
}
