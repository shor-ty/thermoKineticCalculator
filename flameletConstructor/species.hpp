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
#ifndef species_hpp
#define species_hpp
//- system headers
#include <iostream>
//- user def. headers
//#include "constantValues.hpp"
#include "transportProperties.hpp"
#include "thermodynamicProperties.hpp"

class Species : public TransportProperties, public ThermodynamicProperties
{
    public:

        //- constructor
        Species();

        //- destructor
        ~Species();


    public:

    //- functions
    //-----------

        //- set species name
        void setName( const std::string& );

        //- get species name
        std::string getName() const;

        //- set species ID
        void setID( const int& );

        //- get species ID
        int getID() const;

        //- show values
        void showValues() const;

        //- set thermodynamic and transport properties (if matched)
        void setTransportThermodynamicToTrue();

        //- get bool
        bool getBool() const;


    private:

    //- variables
    //-----------

        //- species name (chemical symbols)
        std::string name;

        //- species ID
        int id{0};

        //- thermodynamic and transport properties
        bool thermoAndTransPropAvailable{false};
};


#endif // species_hpp
