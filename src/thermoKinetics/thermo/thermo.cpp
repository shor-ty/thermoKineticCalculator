/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     T hermo       | Open Source Thermo-Kinetic Library
  c-o-o-c     K iknetic     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Constructor.

    TKC is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    TKC is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TKC; if not, see <http://www.gnu.org/licenses/>

\*---------------------------------------------------------------------------*/

#include "thermo.hpp"
#include "constants.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::Thermo::Thermo(const string fileName, const bool thermo)
:
    ThermoCalc(fileName)
{
    if (debug_)
    {
        Info<< "Constructor Thermo\n" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::Thermo::~Thermo()
{
    if (debug_)
    {
        Info<< "Destructor Thermo\n" << endl;
    }
}


// * * * * * * * * * * * * * * * Check Function  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * Summary Function  * * * * * * * * * * * * * * //

void TKC::Thermo::summary(ostream& data) const
{
    //- Header
    data<< Header() << "\n";

    const wordList& tspecies = species();
    const wordList& tformula = formula();

    data<< " c-o Thermodynamic summary:\n"
        << " ==========================\n\n"
        << " Species in thermo: " << tspecies.size() << "\n";

    data<< " \n\n Species used:\n";
    data<< std::left;

    forEach(tspecies, s)
    {
        data<< "    |--> "
            << std::setw(20) <<  tspecies[s]
            << "(" << tformula[s] << ")\n";
    }

    //- Build table with all NASA coeffs for all species
    NASAPolynomials(data, "LOW");
    NASAPolynomials(data, "HIGH");

    //- Build thermoanalyse table
    thermoTable(data);

}


void TKC::Thermo::NASAPolynomials(ostream& data, const word coeff) const
{
    word range;

    if (coeff == "LOW")
    {
        range = "LOW";
    }
    else if (coeff == "HIGH")
    {
        range = "HIGH";
    }

    //- Header
    data<< "\n\n"
        << "=================================================================="
        << "=================================================================="
        << "================\n"
        << " c-o NASA POLYNOMIALS (" << range << " TEMPERATURE)\n"
        << "=================================================================="
        << "=================================================================="
        << "================\n\n"
        << " Species                |"
        << "         c1               c2               c3      "
        << "         c4               c5               c6      "
        << "         c7        |\n"
        << "------------------------------------------------------------------"
        << "------------------------------------------------------------------"
        << "----------------\n";

    //- Species of Thermodynamic Data
    const wordList& tspecies = species();

    //- Build Table
    forEach(tspecies, s)
    {

        List<scalar> NASA(7, 0);

        if (coeff == "LOW")
        {
            NASA = NASACoeffsLT(tspecies[s]);
        }
        else if (coeff == "HIGH")
        {
            NASA = NASACoeffsHT(tspecies[s]);
        }

        //- For species number + name
        std::ostringstream oss;

        oss << " (" << toStr(s+1) << ") " << tspecies[s];

        data<< std::left << std::setw(22) << oss.str()
            << "  |" << std::right
            << std::setw(17) << NASA[0]
            << std::setw(17) << NASA[1]
            << std::setw(17) << NASA[2]
            << std::setw(17) << NASA[3]
            << std::setw(17) << NASA[4]
            << std::setw(17) << NASA[5]
            << std::setw(17) << NASA[6]
            << "  |\n";
    }

    data<< "=================================================================="
        << "=================================================================="
        << "================\n\n";
}


void TKC::Thermo::thermoTable(ostream& data) const
{
    //- Species of Thermodynamic Data
    const wordList& tspecies = species();
    const wordList& tformula = formula();

    //_ Build the thermoanalyse table
    forEach(tspecies, ts)
    {
        const word& tphase = phase(tspecies[ts]);

        data<< "==========================================================="
            << "=========================================================\n"
            << " c-o Thermo analyses for " << tspecies[ts] << "\n"
            << "==========================================================="
            << "=========================================================\n"
            << " Phase: " << tphase << "\n"
            << " Formula: " << tformula[ts] << "\n"
            << " Composition:\n   |\n";

        const map<word, scalar>& elements = elementAtomsMap(tspecies[ts]);

        forAll(elements, a)
        {
            data<< "   |--> " << std::setw(4) << a.first
                << " " << a.second << "\n";
        }

        data<< "\n" << std::setw(40) << std::left
            <<" Molecular weight:   "
            << std::setw(14) << std::right
            << MW(tspecies[ts]) << " [g/mol]\n"
            << std::setw(40) << std::left
            << " Formation enthalpy (298K):   "
            << std::setw(14) << std::right
            << hf(tspecies[ts]) << " [J/mol]\n"
            << std::setw(40) << std::left
            <<" Frormation free Gibbs energy (298K):   "
            << std::setw(14) << std::right
            << gf(tspecies[ts]) << " [J/mol]\n";


        data<< std::right << "\n";

        data<< "\n\n-----------------------------------------------------------"
            << "---------------------------------------------------------\n"
            << "      T    |         cp              H               S       "
            << "        G              dHf            dGf          |\n"
            << "     [K]   |      [J/molK]        [J/mol]         [J/molK]   "
            << "     [J/mol]         [J/molK]       [J/mol]        |\n"
            << "-----------------------------------------------------------"
            << "---------------------------------------------------------\n";

        for(int i=300; i<=3000; i+=100)
        {
            data<< "  " << std::setw(6) << i << "   |"
                << "  " << std::setw(13) << cp(tspecies[ts], i)
                << "  " << std::setw(14) << h(tspecies[ts], i)
                << "  " << std::setw(14) << s(tspecies[ts], i)
                << "  " << std::setw(14) << g(tspecies[ts], i)
                << "  " << std::setw(14) << dhf(tspecies[ts], i)
                << "  " << std::setw(14) << dgf(tspecies[ts], i)
                << "     |\n";
        }

        data<< " ----------------------------------------------------------"
            << "---------------------------------------------------------\n\n"
            << "\n\n";
    }
}


// ************************************************************************* //
