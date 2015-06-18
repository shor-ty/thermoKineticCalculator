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
#include <iomanip>

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
    Chemistry chemistry;

    //- Read kinetic file
    readChemKinThermo(fKinetic, fThermo, species, chemistry);

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
    //scalar Tst_a = adiabateFlameTemperature(Zst, species);

    //- mixture fraction Z (discrete points)
    scalarField Z_dP = discretZ(fAFCDict, Zst);

    //- scalar dissipation rates (discrete points)
    scalarField chi_dP = discretChi(fAFCDict);

    summary(hf_a, ho_a, Zst, chi_dP, species);

    //- Now the complex stuff
    //  for each discrete point Z, store one object species
    std::vector<std::vector<Species> > Z;
    forAll(Z_dP, point)
    {
        Z.push_back(species);
    }

    //- Set up the fuel and oxidizer stream for start (no chemical reaction)
    initializeY(Z_dP, Z);


    std::cout<< "No. \t Reaction \t\t kf   kb \t   LOW   TROE \t    A0  \t    n\t          Ea\n";
    std::cout<< "----------------------------------------------------------------------------------------------------------\n";
    for(int i=0; i<chemistry.r(); i++)
    {
        std::cout<<std::setw(2)  <<  i+1 << "/"
                 << chemistry.r()<<std::setw(3) << ": "
                 << std::left<< std::setw(22) <<chemistry.elementarReaction(i) << " | "
                 << chemistry.fb(i,0) << "    " << chemistry.fb(i,1) << " | \t | "
                 << chemistry.TROE(i) << "      " << chemistry.LOW(i) << " | \t | "
                 << std::setw(9) << chemistry.arrheniusCoeffs(i,0) << "      " << std::right<< std::setw(5) << chemistry.arrheniusCoeffs(i,1) << "      "<< std::setw(8) << chemistry.arrheniusCoeffs(i,2)<< std::setw(5) << " |"
                 << std::endl;

    }

    return 0;
}
