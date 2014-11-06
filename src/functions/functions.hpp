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
#ifndef functions_hpp
#define functions_hpp
//- system headers
//- user def. headers
#include "species.hpp"
#include "chemistry.hpp"

//- function declaration
//----------------------

    //- check input parameters
    void checkInputParameters( const int&, const char*, const char* );

    //- show start info
    void showInfo();

    //- error message; wrong parameters
    void showErrorArgumentInput();

    //- error message; wrong input
    void showErrorWrongInput();

    //- create species objects
    std::vector<Species> createSpeciesObjects( const std::string& );

    //- get transport data
    void getTransportData( const std::string&, std::vector<Species>& );

    //- get thermo data
    void getThermoData( const std::string&, std::vector<Species>& );

    //- get kinetic data
    Chemistry generateAndGetKineticData( const std::string& );

#endif // functions_hpp
