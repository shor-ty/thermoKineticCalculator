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

#include "chemistry.hpp"
#include "thermo.hpp"
#include "constants.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Chemistry::Chemistry(const string fileName, const Thermo& thermo)
:
    thermo_(thermo)
{
    ChemistryReader chemReader(fileName);

    chemReader.read(chemData_);

    //- Reactions that include species i
    createSpeciesInReaction();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Chemistry::~Chemistry()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool AFC::Chemistry::thermo()
{
    return (chemData_.thermo());
}


// * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * * //

AFC::scalar AFC::Chemistry::calculateOmega
(
    const word species,
    const scalar T,
    const map<word, scalar>& con
) const
{
    //- Calculate source term omega
    return chemCalc_.calculateOmega(species, T, con, thermo_, chemData_);
}


AFC::map<AFC::word, AFC::scalar> AFC::Chemistry::calculateOmega
(
    const scalar T,
    const map<word, scalar>& con
) const
{
    //- Species
    const wordList& species = this->species();

    //- Temporary map
    map<word, scalar> dcdt;

    //- Build the map
    forAll(species, s)
    {
        dcdt[s] = scalar(0);
    }

    forAll(species, s)
    {
        //- Calculate source term omega for species s
        dcdt[s] = calculateOmega(s, T, con);
    }

    //- Return the rate field
    return dcdt;
}


AFC::scalar
AFC::Chemistry::kf(const int r, const scalar T, const bool LOW) const
{
    return chemCalc_.kf(r, T, chemData_, LOW);
}


AFC::scalar
AFC::Chemistry::kb(const int r, const scalar T, const bool LOW) const
{
    return chemCalc_.kb(r, T, chemData_, thermo_, LOW);
}


AFC::scalar AFC::Chemistry::keq(const int r, const scalar T) const
{
    return chemCalc_.keq(r, T, chemData_, thermo_);
}


AFC::scalar AFC::Chemistry::Fcent(const int r, const scalar T) const
{
    return chemCalc_.Fcent(r, T, chemData_);
}


AFC::scalar
AFC::Chemistry::Flog(const int r, const scalar T, const scalar M) const
{
    return chemCalc_.Flog(r, T, M, chemData_);
}


void AFC::Chemistry::calculateKb()
{
    //chemCalc_.calculateKb(chemData_);
}


// * * * * * * * * * * * * * * * Update Functions  * * * * * * * * * * * * * //

//void AFC::Chemistry::updateM
//(
//    const scalar& M
//)
//{
//    chemData_.updateM(M);
//}


// * * * * * * * * * * * * * * * Create Functions  * * * * * * * * * * * * * //

void AFC::Chemistry::createSpeciesInReaction()
{
    //- Species list
    const wordList& species = chemData_.species();

    //- Reaction no.
    const int nReac = chemData_.nReac();

    //- List<wordList> that contains all species in each reaction
    const List<wordList>& speciesInReaction = chemData_.speciesInReaction();

    forAll(species, s)
    {
        //- If found species in reaction -> true
        bool found{false};

        //- Loop through all elementar reactions
        for(int r=0; r<nReac; r++)
        {
            //- Loop through the species i in elementar reaction r
            for(unsigned int i=0; i<speciesInReaction[r].size(); i++)
            {
                if
                (
                  ! found
                 && speciesInReaction[r][i] == s
                )
                {
                    found = true;

                    //- Insert reaction no to matrix
                    chemData_.setReacNumbers(s, r);
                }
            }

            //- Reset
            found = false;
        }
    }
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

bool AFC::Chemistry::BR(const int reacNo) const
{
    return chemData_.BR(reacNo);
}


bool AFC::Chemistry::TBR(const int reacNo) const
{
    return chemData_.TBR(reacNo);
}


bool AFC::Chemistry::LOW(const int reacNo) const
{
    return chemData_.LOW(reacNo);
}


bool AFC::Chemistry::TROE(const int reacNo) const
{
    return chemData_.TROE(reacNo);
}


bool AFC::Chemistry::SRI(const int reacNo) const
{
    return chemData_.SRI(reacNo);
}


bool AFC::Chemistry::ENHANCED(const int reacNo) const
{
    return chemData_.ENHANCED(reacNo);
}


AFC::wordList AFC::Chemistry::species() const
{
    return chemData_.species();
}


AFC::wordList AFC::Chemistry::elements() const
{
    return chemData_.elements();
}


unsigned int AFC::Chemistry::nDublicated() const
{
    return chemData_.nDublicated();
}


int AFC::Chemistry::nReac() const
{
    return chemData_.nReac();
}


AFC::string AFC::Chemistry::elementarReaction(const int r) const
{
    return chemData_.elementarReaction(r);
}


AFC::List<AFC::string> AFC::Chemistry::elementarReaction() const
{
    return chemData_.elementarReaction();
}


//AFC::scalarField AFC::Chemistry::reacNoForSpecies
//(
//    const int s
//) const
//{
//    return chemData_.reacNoForSpecies(s);
//}
//
//
//AFC::scalarField AFC::Chemistry::k() const
//{
//    return chemData_.k();
//}*/


//AFC::scalar AFC::Chemistry::k
//(
//    const int reacNo
//) const
//{
//    return chemData_.k(reacNo);
//}


AFC::scalar AFC::Chemistry::dH(const int r, const scalar T) const
{
    return chemCalc_.dH(r, T, chemData_, thermo_);
}


AFC::scalar AFC::Chemistry::dG(const int r, const scalar T) const
{
    return chemCalc_.dG(r, T, chemData_, thermo_);
}


AFC::scalar AFC::Chemistry::dS(const int r, const scalar T) const
{
    return chemCalc_.dS(r, T, chemData_, thermo_);
}


AFC::List<int> AFC::Chemistry::reacNumbers(const word species) const
{
    return chemData_.reacNumbers(species);
}


AFC::wordList AFC::Chemistry::speciesInReaction(const int r) const
{
    return chemData_.speciesInReaction(r);
}


AFC::wordList AFC::Chemistry::speciesProducts(const int r) const
{
    return chemData_.speciesProducts(r);
}


AFC::wordList AFC::Chemistry::speciesEducts(const int r) const
{
    return chemData_.speciesEducts(r);
}


AFC::map<AFC::word, int> AFC::Chemistry::nuProducts(const int r) const
{
    return chemData_.nuProducts(r);
}


AFC::map<AFC::word, int> AFC::Chemistry::nuEducts(const int r) const
{
    return chemData_.nuEducts(r);
}


// * * * * * * * * * * * * * * * Summary Functions * * * * * * * * * * * * * //

void AFC::Chemistry::summary(ostream& data) const
{

    //- Header
    data<< Header() << "\n";

    const wordList& elements = chemData_.elements();
    const wordList& species = chemData_.species();
    const wordList& reactions = chemData_.elementarReaction();

    //- Get amount of LOW TROE SRI ENHANCED BR IR DUB reactions
    unsigned int BR{0};
    unsigned int IR{0};
    unsigned int LOW{0};
    unsigned int TROE{0};
    unsigned int SRI{0};
    unsigned int ENH{0};
    unsigned int TBR{0};

    unsigned int DUB = chemData_.nDublicated();

    forEach(reactions, r)
    {
        if (chemData_.BR(r))
        {
            BR++;
        }
        else
        {
            IR++;
        }
        if (chemData_.LOW(r))
        {
            LOW++;
        }
        if (chemData_.TROE(r))
        {
            TROE++;
        }
        if (chemData_.SRI(r))
        {
            SRI++;
        }
        if (chemData_.ENHANCED(r))
        {
            ENH++;
        }
        if (chemData_.TBR(r))
        {
            TBR++;
        }
    }

    //- General stuff
    data<< " c-o Reaction summary:\n"
        << "======================\n\n"
        << " Number of elements             " << elements.size()  << "\n"
        << " Number of species              " << species.size()   << "\n"
        << " Number of dublicated reaction  " << DUB  << "\n"
        << " Number of elementar reactions  " << reactions.size() << "\n"
        << "   |\n"
        << "   |--> Number of equilibrium reactions   " << BR   << "\n"
        << "   |--> Number of irreversible reactions  " << IR   << "\n"
        << "   |--> Number of ThirdBody reactions     " << TBR  << "\n"
        << "   |--> Number of ENHANCED reactions      " << ENH  << "\n"
        << "   |--> Number of LOW reactions           " << LOW  << "\n"
        << "   |--> Number of TROE reactions          " << TROE << "\n"
        << "   |--> Number of SRI reactions           " << SRI  << "\n";

    data<< "\n\n Elements used: \n   | \n";

    forEach(elements, e)
    {
        data<< "   |-->  (" << e+1 << ")  " << elements[e] << "\n";
    }

    data<< "\n\n Species used:\n   | \n";

    forEach(species, s)
    {
        data<< "   |-->  (" << s+1 << ")  " << species[s] << "\n";
    }

    data<< "\n\n Reactions used: \n   | \n";

    forEach(reactions, r)
    {
        data<< "   |-->  (" << r+1 << ")  " << reactions[r] << "\n";
    }

    data<< "\n\n More detailed kinetics for reactions\n\n";

    chemicalTable(data);

}


void AFC::Chemistry::chemicalTable(ostream& data) const
{
    //- Build the chemical table
    const wordList& reactions = chemData_.elementarReaction();

    List<word> UNITS{ "[1/s]" , "[cm^3/mol/s]", "[cm^6/mol^2/s]" };

    forEach(reactions, r)
    {
        //- Get information of reaction order (forward and backward)
        const scalar& fRO = chemData_.forwardReactionOrder(r);
        const scalar& bRO = chemData_.backwardReactionOrder(r);
        const scalar& globalRO = chemData_.globalReactionOrder(r);

        //- Units of rate constant k (depend on reaction order)
        word unitskf{""};
        word unitskb{""};

        //- TODO maybe a new field in chemData
        if (fRO == scalar(-1))
        {
            unitskf = UNITS[0];
        }
        else if (fRO == scalar(-2))
        {
            unitskf = UNITS[1];
        }
        else if (fRO == scalar(-3))
        {
            unitskf = UNITS[2];
        }

        if (bRO == scalar(1))
        {
            unitskb = UNITS[0];
        }
        else if (bRO == scalar(2))
        {
            unitskb = UNITS[1];
        }
        else if (bRO == scalar(3))
        {
            unitskb = UNITS[2];
        }

        const List<scalar>& arrheniusCoeffs = chemData_.arrheniusCoeffs(r);

        data<< "========================================================"
            << "=====================================================\n"
            << " c-o  Reaction " << r+1 << ":  " << reactions[r] << " at "
            << thermo_.p() << " Pa\n"
            << "========================================================"
            << "=====================================================\n"
            << std::left  << std::setw(30) << "   Reaction order forward:"
            << std::right << std::setw(15) << fRO << "\n"
            << std::left  << std::setw(30) << "   Reaction order backward:"
            << std::right << std::setw(15) << bRO << "\n"
            << std::left  << std::setw(30) << "   Change in moles:"
            << std::right << std::setw(15) << globalRO << "\n"
            << std::left  << std::setw(30) << "   Units of k (forward): "
            << std::right << std::setw(15) << unitskf << "\n"
            << std::left  << std::setw(30) << "   Units of k (backward): "
            << std::right << std::setw(15) << unitskb << "\n"
            << "\n\n"
            << "   Arrhenius coefficients\n"
            << "--------------------------------------------------\n   |\n"
            << "   |-->  A:  " << std::setw(14) << arrheniusCoeffs[0] << " "
            << "\n"
            << "   |-->  n:  " << std::setw(14) << arrheniusCoeffs[1] << " "
            << "\n"
            << "   |-->  Ea: " << std::setw(14) << arrheniusCoeffs[2] << " "
            << "cal/mol\n\n";

        if(chemData_.TBR(r))
        {

            if(chemData_.LOW(r))
            {
                data<< "   LOW Coeffs (for high pressure)\n"
                    << "--------------------------------------------------\n"
                    << "   |\n";

                const List<scalar>& LOWCoeffs = chemData_.LOWCoeffs(r);

                data<< "   |-->  A:  " << std::setw(14) << LOWCoeffs[0] << " ";
                data<< "\n";
                data<< "   |-->  n:  " << std::setw(14) << LOWCoeffs[1] << " ";
                data<< "\n";
                data<< "   |-->  Ea: " << std::setw(14) << LOWCoeffs[2] << " ";
                data<< "cal/mol\n\n";
            }

            if(chemData_.TROE(r))
            {
                data<< "   TROE Coeffs\n"
                    << "--------------------------------------------------\n"
                    << "   |\n";

                const List<scalar>& TROECoeffs = chemData_.TROECoeffs(r);

                data<< "   |-->  a:  " << TROECoeffs[0] << "\n";
                data<< "   |-->  b:  " << TROECoeffs[1] << "\n";
                data<< "   |-->  c:  " << TROECoeffs[2] << "\n";
                data<< "   |-->  d:  " << TROECoeffs[3] << "\n\n";
            }

            if(chemData_.ENHANCED(r))
            {
                data<< "   Third-Body Coeffs adjustment\n"
                    << "--------------------------------------------------\n"
                    << "   |\n";

                const map<word, scalar>& enhanced = chemData_.ENHANCEDCoeffs(r);

                forAll(enhanced, species)
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

        if (chemData_.LOW(r))
        {
            buildTablekf(r, data, true);
        }

        if (chemData_.TROE(r))
        {
            buildTROETable(r, data);
        }
    }
}


void AFC::Chemistry::buildTablekf
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
            << "  " << std::setw(13) << dH(r, i)
            << "  " << std::setw(13) << dG(r, i)
            << "  " << std::setw(13) << dS(r, i) << "   |\n";
    }

    data<< "--------------------------------------------------------"
        << "-----------------------------------------------------\n"
        << "\n\n";
}


void AFC::Chemistry::buildTROETable(const int r, ostream& data) const
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


// ************************************************************************* //
