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

//#include "thermodynamic/thermodynamic.hpp"
#include "functions/functions.hpp"
#include "database/elements.hpp"


int main()
{
    normalString s = "../files/kinetics.kin";
    normalString t = "../files/thermo.tdc";
    //- Some output definitions
//    std::cout.precision( 8 );
    //std::cout.setf( std::ios::scientific );


    //- Read kinetic file
    readChemKinThermo(s, t);




    //- read thermodynamic file
//    const stringField thermoFileContent = openFile("../files/thermodynamic.tdc");

    //- create thermodynamic objects
//    std::vector<Thermodynamic> thermo = createThermodynamicObjects(thermoFileContent);

    //for(int i=0; i< 1000; i=i+50)

//    Info << "buba" ;
    return 0;
}
