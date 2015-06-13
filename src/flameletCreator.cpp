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
    normalString fKinetic = "../files/kinetics.kin";
    normalString fThermo  = "../files/thermo.tdc";
    normalString fAFCDict = "../files/afcDict";
    //- Some output definitions
    //std::cout.precision( 8 );
    //std::cout.setf( std::ios::scientific );

    //- SPECIES class
    std::vector<Species> species;

    //- REACTION class
    std::vector<Reactions> reactions;

    //- Read kinetic file
    readChemKinThermo(fKinetic, fThermo, species, reactions);

    //- Read afcDict
    readAFCDict(fAFCDict, species);

    //- calculate stochiometric mixture fraction Zst
    scalar Zst = stochiometricMF(species);

    //- adiabatic enthalpy of fuel and oxidizer
    //  + 0 mean fuel
    //  + 1 mean oxidizer
    //  [J/kg]
    scalar hf_a = adiabaticEnthalpy(species, 0);
    scalar ho_a = adiabaticEnthalpy(species, 1);

    //- calculate adiabatic flame temperature
    //  for stochiometric conditions
    scalar Tst_a = adiabateFlameTemperature(Zst, species);

    //- mixture fraction Z (discrete points)
    scalarField Z_dP = discretZ(fAFCDict);

    //- scalar dissipation rates (discrete points)
    scalarField chi_dP = discretChi(fAFCDict);

    summary(hf_a, ho_a, Zst, chi_dP, species);


    //- partial differential equation
    //scalar t{0};

    return 0;
}
