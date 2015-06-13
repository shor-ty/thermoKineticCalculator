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

//- set mol fraction X in fuel
void Species::setXf(const scalar& molFraction)
{
    Xf_ = molFraction;
}

//- set mol fraction X in oxidizer
void Species::setXo(const scalar& molFraction)
{
    Xo_ = molFraction;
}

//- return mol fraction X in fuel
const scalar Species::Xf() const
{
    return Xf_;
}

//- return mol fraction X in oxidizer
const scalar Species::Xo() const
{
    return Xo_;
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
