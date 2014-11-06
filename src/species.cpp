/*---------------------------------------------------------------------------*\
| =========                |                                                  |
| \\      /  A utomatic    | Holzmann-cfd                                     |
|  \\    /   F lamelet     | Version: 1.0                                     |
|   \\  /    C constructor | Web: www.Holzmann-cfd.de                         |
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
#include "species.hpp"

Species::Species()
{
 std::cout << "Object created\n";
}

Species::~Species()
{
 std::cout << "Object destroyed\n";
}

//- functions
//-----------

    void Species::setName( const std::string& nameOfSpecies )
    {
        name = nameOfSpecies;
    }

    void Species::setID( const int& speciesID )
    {
        id = speciesID;
    }

    void Species::setTransportThermodynamicToTrue()
    {
        thermoAndTransPropAvailable = true;
    }

    void Species::showValues() const
    {
        std::cout << name << "\n\n";

        //-
        TransportProperties::showValues();

        //-
        ThermodynamicProperties::showValues();

        std::cout << "------------------------------------------------------------\n";
        std::cout << std::endl;
    }

//- functions
//-----------

    std::string Species::getName() const
    {
        return name;
    }

    int Species::getID() const
    {
        return id;
    }

    bool Species::getBool() const
    {
        return thermoAndTransPropAvailable;
    }

