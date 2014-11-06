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
#include <iostream>
#include <vector>
#include <string>
#include <memory>

//- header files
#include "functions.hpp"
#include "species.hpp"
#include "chemistry.hpp"
//#include "constantValues.hpp"


//#include <boost/algorithm/string.hpp>
//#include <math.h>
//#include <stdlib.h>
//#include <limits>


int main( int argc, char* argv[] )
{
    //- Some output definitions
    std::cout.precision( 9 );
    std::cout.setf( std::ios::scientific );

    //-
    showInfo();

    //- file names
    std::string filePath_kin, filePath_tdc, filePath_tra;

    //- check input parameters
    for( int i=1; i<argc; i+=2 )
    {
        checkInputParameters( argc, argv[i], argv[i+1] );
    }

    //- get file name
    for( int i=1; i<argc; i+=2 )
    {
        std::string arg_1 = argv[i];

        if( arg_1 == "-thermo" )     { std::cout << " + tdc file \t" << argv[i+1] << std::endl; filePath_tdc = std::string(argv[i+1]); }
        if( arg_1 == "-kinetic" )    { std::cout << " + kin file \t" << argv[i+1] << std::endl; filePath_kin = std::string(argv[i+1]); }
        if( arg_1 == "-transport" )  { std::cout << " + tra file \t" << argv[i+1] << std::endl; filePath_tra = std::string(argv[i+1]); }
    }

    //- generate an vector array of all species objects
    auto vecOfAllSpecies = createSpeciesObjects( filePath_tra );

    //- get transport data
    getTransportData( filePath_tra, vecOfAllSpecies );

    //- get thermodynamic data
    getThermoData( filePath_tdc, vecOfAllSpecies );

    //- generate the chemistry object and get kinetic data
    auto chemistryObj = generateAndGetKineticData( filePath_kin );


    //chemistryObj.showSummary();

    std::cout << "\n + Programm successfully closed...\n\n";
    return 0;
}

