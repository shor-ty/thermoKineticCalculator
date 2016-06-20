/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2015 Holzmann-cfd
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Constructor.

    AFC is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or 
    (at your option) any later version.

    AFC is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AFC; if not, see <http://www.gnu.org/licenses/>

\*---------------------------------------------------------------------------*/

#include <fstream>
#include <iomanip>
#include "typedef.hpp"
#include "constants.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

class Thermo;
class Chemistry;

// * * * * * * * * * * * * * AFC Interprete function * * * * * * * * * * * * //

void interprete
(
    Thermo thermo,
    Chemistry chemistry,   
    Properties properties
)
{
    //- Create folder TODO - use BOOST for cross platforms
    system("mkdir -p interpretedData");

    //- Interprete the thermo data
    wordList species;
    
        //- Which species should be analysed in a thermodynamic way?
        if (properties.interpreter() == "THERMO")
        {
            Info<< "    c-o Analyse all species in thermo database \n";
            species = thermo.species();
        }
        else
        {
            Info<< "    c-o Analyse all species in chemistry database \n";
            species = chemistry.species();
        }

    //- File to write in
    std::ofstream output;

    output.open
    (   
        "interpretedData/thermoAnalyse",
        std::ofstream::out | std::ofstream::trunc
    );

    Info<< "    c-o Thermodynamic analyses ...";

    forAll(species, s)
    {

        const word& speci = species[s];

        output << " *** THERMO ANALYSE OF " + speci + "\n"
            << "---------------------------------------------------------\n"
            << "  c-o Phase: " + thermo.phase(speci) + "\n"
            << "  c-o Composition: ";

        //- Composition
        const wordList& atoms = thermo.elementsInSpecies(speci);
        const scalarList& factors = thermo.elementsFactors(speci);

        forAll(atoms, a)
        {
            output << factors[a] << atoms[a] << "  ";
        }



        output << "\n  c-o Moleculare weight: " << thermo.MW(speci) << " g/mol\n"
            << "  c-o Free Gibbs energy at 298 K: " <<  thermo.G(speci, 298)
            << " J/mol\n"
            << "  c-o Formation enthalpy at 298 K: " <<  thermo.Hf(speci, 298)
            << " J/mol\n"
            << "==================================================="
            << "===========================================\n"
            << std::setw(12) << "Temperature" 
            << std::setw(8) << "Cp"
            << std::setw(12) << "H"
            << std::setw(12) << "S"
            << std::setw(12) << "G"
            << std::setw(12) << "dH"
            << std::setw(12) << "dG"
            << "\n    "
            << std::setw(4) << "[K]"
            << "    "
            << std::setw(12) << "[J/mol/K]"
            << std::setw(11) << "[J/mol]"
            << std::setw(13) << "[J/mol/K]"
            << std::setw(12) << "[J/mol/K]"
            << std::setw(10) << "[J/K]"
            << std::setw(11) << "[J/K]"
            << "\n"
            << "==================================================="
            << "===========================================\n";

        const scalar& TL = thermo.LT(speci);
        const scalar& TH = thermo.HT(speci);

        for (scalar T = TL; T <= TH;)
        {
            output << "    " << std::setw(12) << std::left <<  T 
                << std::setw(12) << std::left << thermo.cp(speci, T)
                << std::setw(12) << std::left << thermo.H(speci, T)
                << std::setw(12) << std::left << thermo.S(speci, T) 
                << std::setw(12) << std::left << thermo.G(speci, T) 
                << std::setw(12) << std::left << thermo.H(speci, T) - thermo.H(speci, 298)
                << std::setw(12) << std::left << thermo.G(speci, T) - thermo.G(speci, 298)
                << "  \n";
            T += 100;
        }

        output << "==================================================="
            << "===========================================\n\n\n";
    }

    Info<< " done\n";

    //- Close file for thermo analyse
    output.close();


    //- Analyse the chemical stuff
    Info<< "    c-o Chemistry analyse ...";
    
    //- Interprete the thermo data
    const wordList& reactions = chemistry.elementarReaction();

    //- File to write in
    output.open
    (   
        "interpretedData/chemistryAnalyse",
        std::ofstream::out | std::ofstream::trunc
    );

    forAll(reactions, r)
    {

        const bool& TBR = chemistry.TBR(r);
        const bool& ENHANCED = chemistry.ENHANCED(r);
        const bool& LOW = chemistry.LOW(r);
        const bool& TROE = chemistry.TROE(r);
        const bool& SRI = chemistry.SRI(r);

        //- Arrhenius coeffs
        const scalarList& arrheniusCoeffs = chemistry.arrheniusCoeffs(r);
        const scalar& A = arrheniusCoeffs[0];
        const scalar& beta = arrheniusCoeffs[1];
        const scalar& Ea = arrheniusCoeffs[2];

        output << " *** ELEMENTAR REACTION " << r+1 << ": " << reactions[r] << "\n"
            << "----------------------------------------------------------------\n";

        output << "\n  c-o Pressure: " << thermo.p() << " [Pa]";

        if (TBR)
        {

            output<< "\n";

            if (!TROE && !LOW && !SRI)
            {
                output <<std::setw(13)<< "  c-o kf: " << A << " * T^" << beta
                    << " * e^(" << -Ea << "/(RT))" << "   [cm^3/mol/s]\n";
            }

            if (LOW)
            {
                output << "  c-o Fall-off reaction\n";
                output << "  c-o Thirdbody [M] == [" << properties.inertGas() << "]\n"
                       << "  -----------------------------------\n";

                //- Arrhenius coeffs LOW pressure
                const scalarList& arrheniusCoeffs = chemistry.LOWCoeffs(r);
                const scalar& ALow = arrheniusCoeffs[0];
                const scalar& betaLow = arrheniusCoeffs[1]; 
                const scalar& EaLow = arrheniusCoeffs[2];

                //- Low 
                output <<std::setw(17)<< "  c-o kf (low): " << ALow << " * T^" << betaLow
                    << " * e^(" << -EaLow << "/(RT))" << "   [cm^3/mol/s]\n";
                //- High
                output <<std::setw(16)<< "  c-o kf (high): " << A << " * T^" << beta
                    << " * e^(" << -Ea << "/(RT))" << "   [cm^3/mol/s]\n";
            }
            if (TROE)
            {
                output << "\n  c-o TROE coeffs:\n"
                       << "  -----------------------------------\n";

                //- TROE coeffs
                const scalarList& TROECoeffs = chemistry.TROECoeffs(r);
                const scalar& alpha = TROECoeffs[0];
                const scalar& Tsss = TROECoeffs[1];
                const scalar& Ts = TROECoeffs[2];
                const scalar& Tss = TROECoeffs[3];

                output << "      a     = " << std::setw(10) << alpha << "\n"
                    << "      T***  = " << std::setw(10) << Tsss << "[K]\n" 
                    << "      T*    = " << std::setw(10) << Ts << "[K]\n" 
                    << "      T**   = " << std::setw(10) << Tss << "[K]\n";

            }
            if (ENHANCED)
            {
                output << "\n  c-o Three-Body enhanced factors:\n"
                       << "  -----------------------------------\n";

                const wordList& species = chemistry.enhancedSpecies(r);
                const map<word, scalar> enhancedFactors = chemistry.enhancedFactors(r);

                forAll(species, s)
                {
                    output << "      " << std::setw(5) << species[s]
                        << std::setw(5) 
                        << enhancedFactors.at(species[s]) << "\n";
                }
            }

        }
        else
        {
            output <<std::setw(13)<< "\n  c-o kf: " << A << " * T^" << beta
                << " * e^(" << -Ea << "/(RT))" << "   [cm^3/mol/s]\n";
        }
        output << "==================================================="
            << "===================================================\n"
            << "Temperature       kf               kc                kb     "
            << "       dH           dS          dG     \n"
            << "    [K]       [mol/cm^3/s]         [-]          [mol/cm^3/s]  "
            << "   [J/mol]     [J/mol/K]    [J/mol] \n" 
            << "==================================================="
            << "===================================================\n";

        for (scalar T = 300.; T < 4000.;)
        {

            if (TBR)
            {
                //- For TBR we need the Third-Body Collision concentration [M] 
                //  Y*T*MW [K*g/mol], only inert gas are used 
                scalar YTMW = 1. * T / thermo.MW(properties.inertGas());

                //- Pressure [Pa]
                const scalar& p = properties.p();

                //- Store inert gas concentration as [M] 
                //  Change unit to [mol/cm^3] factor 100cm * 100cm * 100cm / m^3
                chemistry.updateM
                ( 
                    (((p * 1. ) / thermo.MW(properties.inertGas()))
                        / (AFC::Constants::R * YTMW) / 1e6)
                );
            }

            chemistry.calculateKf(r, T);
            chemistry.calculateKc(r, T, thermo);
            chemistry.calculateKb();

            const scalar& kf = chemistry.kf();             
            const scalar& kb = chemistry.kb();
            const scalar& Kc = chemistry.Kc();
            const scalar& dH = chemistry.dH()* AFC::Constants::jouleToCal / 1000;
            const scalar& dS = chemistry.dS()* AFC::Constants::jouleToCal ;
            const scalar& dG = chemistry.dG()* AFC::Constants::jouleToCal / 1000;

            std::setprecision(2);
            output << "    " << std::setw(9) << std::left <<  T 
          <<  "  " << std::setw(15) << kf  << "  " << std::setw(15) << Kc
                << "  " << std::setw(15) << kb  << "  " << std::setw(10) << dH 
                << "  " << std::setw(10) << dS  << "  " << std::setw(10) << dG 
                << "  \n";
            T += 100;
        }

        output << "==================================================="
            << "===================================================\n\n\n";

    }

    Info<< " done\n";

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
