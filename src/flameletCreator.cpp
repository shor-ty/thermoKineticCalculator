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
#include <iostream>
#include <fstream>
//- user def. headers
#include "thermodynamic/thermodynamic.hpp"
#include "functions/functions.hpp"

int main()
{
    //- Some output definitions
    std::cout.precision( 8 );
    //std::cout.setf( std::ios::scientific );

    //- runtime include files
    #include "definitions/IOStream.hpp"

    //- read thermodynamic file
    const stringField thermoFileContent = openFile("../files/thermodynamic.tdc");

    //- create thermodynamic objects
    std::vector<Thermodynamic> thermo = createThermodynamicObjects(thermoFileContent);

    //for(int i=0; i< 1000; i=i+50)
    int i = 298;
//    std::cout << "H bei " << i << "K" << " = " << (thermo[0].calculateEnthalpy(i))/1000 << "\n";
//    std::cout << "C bei " << i << "K" << " = " << (thermo[0].calculateHeatCapacity(i)) << "\n";
//    std::cout << "S bei " << i << "K" << " = " << (thermo[0].calculateEntropy(i)) << "\n";

//    Info << "buba" ;
    return 0;
}
