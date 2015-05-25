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
Species::Species() : MW_(0), Y_(0), X_(0), Con_(0)
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

const scalar Species::MW()
{
    return MW_;
}

//- return name
const normalString Species::name()
{
    return name_;
}
