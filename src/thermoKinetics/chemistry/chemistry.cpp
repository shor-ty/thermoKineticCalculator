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

#include "chemistry.hpp"
#include "thermo.hpp"
#include "constants.hpp"
#include <algorithm>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::Chemistry::Chemistry(const string fileName, const Thermo& thermo)
:
    ChemistryCalc(fileName, thermo)
{
    //- Check if all species are available
    checkSpecies();

    //- Build the table that contains in which reaction each species is included
    buildSpeciesInReactionTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::Chemistry::~Chemistry()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void TKC::Chemistry::checkSpecies() const
{
    Info<< " c-o Checking if chemistry species are available in the NASA"
        << " database..." << endl;

    const wordList& thermoSpecies = thermo().species();
    const wordList& chemistrySpecies = species();
    wordList::const_iterator it;

    forAll(chemistrySpecies, cSpecies)
    {

        it = find(thermoSpecies.begin(), thermoSpecies.end(), cSpecies);

        if (it == thermoSpecies.end())
        {
            ErrorMsg
            (
                "    Species '" + cSpecies + "' of the chemistry is not "
                " available in the thermodynamic database (NASA Polynomials)",
                __FILE__,
                __LINE__
            );
        }
    }

    //- If passed everything is fine
    Info<< "      >> Everything is fine. Proceed...\n" << endl;
}


// * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * * //

/*
TKC::map<TKC::word, TKC::scalar> TKC::Chemistry::omega
(
    const scalar T,
    const map<word, scalar>& con
) const
{
    //- Species
    const wordList& tspecies = species();

    //- Temporary map
    map<word, scalar> dcdt;

    //- Build the map
    forAll(tspecies, s)
    {
        dcdt[s] = scalar(0);
    }

    forAll(tspecies, s)
    {
        //- Calculate source term omega for species s
        dcdt[s] = ChemistryCalc::omega(s, T, con);
    }

    //- Return the rate field
    return dcdt;
}
*/


// * * * * * * * * * * * * * * * Update Functions  * * * * * * * * * * * * * //

//void TKC::Chemistry::updateM
//(
//    const scalar& M
//)
//{
//    updateM(M);
//}


// * * * * * * * * * * * * * * * Create Functions  * * * * * * * * * * * * * //

void TKC::Chemistry::buildSpeciesInReactionTable()
{
    //- Species list
    const wordList& tspecies = species();

    //- Reaction no.
    const int tnReac = nReac();

    //- List<wordList> that contains all species in each reaction
    const List<wordList>& tspeciesInReaction = speciesInReaction();

    forAll(tspecies, s)
    {
        //- If found species in reaction -> true
        bool found{false};

        //- Loop through all elementar reactions
        for(int r=0; r<tnReac; r++)
        {
            //- Loop through the species i in elementar reaction r
            for(unsigned int i=0; i<tspeciesInReaction[r].size(); i++)
            {
                if (!found && (tspeciesInReaction[r][i] == s))
                {
                    found = true;

                    //- Insert reaction no to matrix
                    setReacNumbers(s, r);
                }
            }

            //- Reset
            found = false;
        }
    }
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * Summary Functions * * * * * * * * * * * * * //

/*
void TKC::Chemistry::summary(ostream& data) const
{

    //- Header
    data<< Header() << "\n";

    const wordList& telements = elements();
    const wordList& tspecies = species();
    const wordList& treactions = elementarReaction();
    const wordList& tignoredReactions = ignoredElementarReaction();

    //- Get amount of LOW TROE SRI ENHANCED BR IR DUB reactions
    unsigned int tBR{0};
    unsigned int tIR{0};
    unsigned int tLOW{0};
    unsigned int tTROE{0};
    unsigned int tSRI{0};
    unsigned int tENH{0};
    unsigned int tTBR{0};

    unsigned int tDUB = nDublicated();
    unsigned int tIGN = nIgnored();

    forEach(treactions, r)
    {
        if (BR(r))
        {
            tBR++;
        }
        else
        {
            tIR++;
        }
        if (LOW(r))
        {
            tLOW++;
        }
        if (TROE(r))
        {
            tTROE++;
        }
        if (SRI(r))
        {
            tSRI++;
        }
        if (ENHANCED(r))
        {
            tENH++;
        }
        if (TBR(r))
        {
            tTBR++;
        }
    }

    //- General stuff
    data<< " c-o Reaction summary:\n"
        << "======================\n\n"
        << " Number of elements             " << telements.size() << "\n"
        << " Number of species              " << tspecies.size() << "\n"
        << " Number of dublicated reaction  " << tDUB  << "\n"
        << " Number of ignored reaction     " << tIGN  << "\n"
        << " Number of elementar reactions  " << treactions.size() << "\n"
        << "   |\n"
        << "   |--> Number of equilibrium reactions   " << tBR << "\n"
        << "   |--> Number of irreversible reactions  " << tIR << "\n"
        << "   |--> Number of ThirdBody reactions     " << tTBR << "\n"
        << "   |--> Number of ENHANCED reactions      " << tENH << "\n"
        << "   |--> Number of LOW reactions           " << tLOW << "\n"
        << "   |--> Number of TROE reactions          " << tTROE << "\n"
        << "   |--> Number of SRI reactions           " << tSRI << "\n";

    data<< "\n\n Elements used: \n   | \n";

    forEach(telements, e)
    {
        data<< "   |-->  (" << e+1 << ")  " << telements[e] << "\n";
    }

    data<< "\n\n Species used:\n   | \n";

    forEach(species, s)
    {
        data<< "   |-->  (" << s+1 << ")  " << tspecies[s] << "\n";
    }

    data<< "\n\n Reactions used: \n   | \n";

    forEach(reactions, r)
    {
        data<< "   |-->  (" << r+1 << ")  " << treactions[r] << "\n";
    }

    data<< "\n\n Reactions ignored: \n   | \n";

    forEach(ignoredReactions, r)
    {
        data<< "   |-->  (" << r+1 << ")  " << tignoredReactions[r] << "\n";
    }

    data<< "\n\n More detailed kinetics for reactions\n\n";

    chemicalTable(data);

}


void TKC::Chemistry::chemicalTable(ostream& data) const
{
    //- Build the chemical table
    const wordList& treactions = elementarReaction();

    List<word> UNITS{ "[1/s]" , "[cm^3/mol/s]", "[cm^6/mol^2/s]" };

    forEach(treactions, r)
    {
        //- Get information of reaction order (forward and backward)
        const scalar& tfRO = forwardReactionOrder(r);
        const scalar& tbRO = backwardReactionOrder(r);
        const scalar& tglobalRO = globalReactionOrder(r);

        //- Units of rate constant k (depend on reaction order)
        word unitskf{""};
        word unitskb{""};

        //- TODO maybe a new field in chemData
        if (tfRO == scalar(-1))
        {
            unitskf = UNITS[0];
        }
        else if (tfRO == scalar(-2))
        {
            unitskf = UNITS[1];
        }
        else if (tfRO == scalar(-3))
        {
            unitskf = UNITS[2];
        }

        if (tbRO == scalar(1))
        {
            unitskb = UNITS[0];
        }
        else if (tbRO == scalar(2))
        {
            unitskb = UNITS[1];
        }
        else if (tbRO == scalar(3))
        {
            unitskb = UNITS[2];
        }

        const List<scalar>& tarrCoeffs = arrheniusCoeffs(r);

        data<< "========================================================"
            << "=====================================================\n"
            << " c-o  Reaction " << r+1 << ":  " << treactions[r] << " at "
            << thermo().p() << " Pa\n"
            << "========================================================"
            << "=====================================================\n"
            << std::left  << std::setw(30) << "   Reaction order forward:"
            << std::right << std::setw(15) << tfRO << "\n"
            << std::left  << std::setw(30) << "   Reaction order backward:"
            << std::right << std::setw(15) << tbRO << "\n"
            << std::left  << std::setw(30) << "   Change in moles:"
            << std::right << std::setw(15) << tglobalRO << "\n"
            << std::left  << std::setw(30) << "   Units of k (forward): "
            << std::right << std::setw(15) << unitskf << "\n"
            << std::left  << std::setw(30) << "   Units of k (backward): "
            << std::right << std::setw(15) << unitskb << "\n"
            << "\n\n"
            << "   Arrhenius coefficients\n"
            << "--------------------------------------------------\n   |\n"
            << "   |-->  A:  " << std::setw(14) << tarrCoeffs[0] << " "
            << "\n"
            << "   |-->  n:  " << std::setw(14) << tarrCoeffs[1] << " "
            << "\n"
            << "   |-->  Ea: " << std::setw(14) << tarrCoeffs[2] << " "
            << "cal/mol\n\n";

        if(TBR(r))
        {

            if(LOW(r))
            {
                data<< "   LOW Coeffs (for high pressure)\n"
                    << "--------------------------------------------------\n"
                    << "   |\n";

                const List<scalar>& tLOWCoeffs = LOWCoeffs(r);

                data<< "   |-->  A:  " << std::setw(14) << tLOWCoeffs[0] << " ";
                data<< "\n";
                data<< "   |-->  n:  " << std::setw(14) << tLOWCoeffs[1] << " ";
                data<< "\n";
                data<< "   |-->  Ea: " << std::setw(14) << tLOWCoeffs[2] << " ";
                data<< "cal/mol\n\n";
            }

            if(TROE(r))
            {
                data<< "   TROE Coeffs\n"
                    << "--------------------------------------------------\n"
                    << "   |\n";

                const List<scalar>& tTROECoeffs = TROECoeffs(r);

                data<< "   |-->  a:  " << tTROECoeffs[0] << "\n";
                data<< "   |-->  b:  " << tTROECoeffs[1] << "\n";
                data<< "   |-->  c:  " << tTROECoeffs[2] << "\n";
                data<< "   |-->  d:  " << tTROECoeffs[3] << "\n\n";
            }

            if(ENHANCED(r))
            {
                data<< "   Third-Body Coeffs adjustment\n"
                    << "--------------------------------------------------\n"
                    << "   |\n";

                const map<word, scalar>& tenhanced = ENHANCEDCoeffs(r);

                forAll(tenhanced, species)
                {
                    data<< "   |-->  " << std::left << std::setw(12)
                        << species.first
                        << "   " << species.second << "\n";
                }
                data<< "\n";
            }
        }

        data<< std::right ;

        buildTablekf(r, data);

        if (LOW(r))
        {
            buildTablekf(r, data, true);
        }

        if (TROE(r))
        {
            buildTROETable(r, data);
        }
    }
}


void TKC::Chemistry::buildTablekf
(
    const int r,
    ostream& data,
    const bool LOW
) const
{
    if (LOW)
    {
        data<< "  Table for LOW pressure coefficients\n";
    }

    data<< "--------------------------------------------------------"
        << "-----------------------------------------------------\n"
        << "      T    |"
        << "        kf             kb             Keq   "
        << "         dH            dG              dS        |\n"
        << "     [K]   |"
        << "       [s.a]          [s.a]           [-]   "
        << "      [J/molK]       [J/mol]        [J/molK]     |\n"
        << "--------------------------------------------------------"
        << "-----------------------------------------------------\n"
        << "";

    for (int i=300; i<=3000; i+=100)
    {
        data<< "  " << std::setw(6) << i << "   |"
            << "  " << std::setw(13) << kf(r, i, LOW) << ""
            << "  " << std::setw(13) << kb(r, i, LOW)
            << "  " << std::setw(13) << keq(r, i)
            << "  " << std::setw(13) << dh(r, i)
            << "  " << std::setw(13) << dg(r, i)
            << "  " << std::setw(13) << ds(r, i) << "   |\n";
    }

    data<< "--------------------------------------------------------"
        << "-----------------------------------------------------\n"
        << "\n\n";
}


void TKC::Chemistry::buildTROETable(const int r, ostream& data) const
{
    data<< "  Table for TROE, [M] assumed to be 1 \n"
        << "--------------------------------------------------\n"
        << "      T    |       Fcent           logF       |\n"
        << "     [K]   |        [-]            [-]        |\n"
        << "--------------------------------------------------\n";

    for (int i=300; i<=3000; i+=100)
    {
        data<< "  " << std::setw(6) << i << "   |"
            << "  " << std::setw(13) << Fcent(r, i)
            << "  " << std::setw(13) << Flog(r, i, scalar(1))
            << "    |\n";
    }

    data<< "--------------------------------------------------\n\n";


}
*/

// ************************************************************************* //
