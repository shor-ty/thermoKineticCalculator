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

#include "thermo.hpp"
#include "constants.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Thermo::Thermo
(
    const string& fileName,
    const bool& thermo
)
:
    thermoData_(thermo)
{
    ThermoReader thermoReader(fileName);

    thermoReader.read(thermoData_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Thermo::~Thermo()
{}


// * * * * * * * * * * * * * * * Insert Functions  * * * * * * * * * * * * * //

void AFC::Thermo::p
(
    const scalar& pressure
)
{
    thermoData_.p(pressure);
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

AFC::wordList AFC::Thermo::species() const
{
    return thermoData_.species();
}


AFC::scalar AFC::Thermo::MW
(
    const word& species
) const
{
    return thermoData_.MW(species);
}


AFC::map<AFC::word, AFC::scalar> AFC::Thermo::MW() const
{
    return thermoData_.MW();
}


AFC::scalar AFC::Thermo::MmeanX
(
    const map<word, scalar>& X
) const
{
    return thermoCalc_.MmeanX(X, MW());    
}


AFC::scalar AFC::Thermo::cp
(
    const word& species,
    const scalar& T
) const
{
    return thermoCalc_.cp(species, T, thermoData_);
}


AFC::scalar AFC::Thermo::H
(
    const word& species,
    const scalar& T
) const
{
    return thermoCalc_.H(species, T, thermoData_);
}


AFC::scalar AFC::Thermo::S
(
    const word& species,
    const scalar& T
) const
{
    return thermoCalc_.S(species, T, thermoData_);
}


AFC::scalar AFC::Thermo::G
(
    const word& species,
    const scalar& T
) const
{
    return thermoCalc_.G(species, T, thermoData_);
}


AFC::scalar AFC::Thermo::G
(
    const scalar& H,
    const scalar& S,
    const scalar& T
) const
{
    return thermoCalc_.G(H, S, T);
}


AFC::scalar AFC::Thermo::p() const
{
    return thermoData_.p();
}


AFC::scalar AFC::Thermo::Hf
(
    const word& species
) const
{
    return thermoCalc_.Hf(species, thermoData_);
}


AFC::scalar AFC::Thermo::Gf
(
    const word& species
) const
{
    return thermoCalc_.Gf(species, thermoData_);
}


AFC::scalar AFC::Thermo::dH
(
    const word& species,
    const scalar& T
) const
{
    return thermoCalc_.dH(species, T, thermoData_);
}


AFC::scalar AFC::Thermo::dG
(
    const word& species,
    const scalar& T
) const
{
    return thermoCalc_.dG(species, T, thermoData_);
}


AFC::scalar AFC::Thermo::C
(
    const scalar& T
) const
{
    return thermoCalc_.C(p(), T);
}


AFC::word AFC::Thermo::phase
(
    const word& species
) const
{
    return thermoData_.phase(species);
}

// * * * * * * * * * * * * * * Summary Function  * * * * * * * * * * * * * * //


void AFC::Thermo::summary
(
    ostream& data
) const
{
    const wordList& species = thermoData_.species();
    const wordList& formula = thermoData_.formula();

    const scalar& M = AFC::Constants::jouleToCal;

    data<< " c-o Thermodynamic summary:\n"
        << " ==========================\n\n"
        << " Species in thermo: " << species.size() << "\n";

    data<< " \n\n Species used:\n";
    data<< std::left;

    forEach(species, s)
    {
        data<< "    |--> " << std::setw(20) <<  species[s]
            << "(" << formula[s] << ")\n";
    }

    forEach(species, s)
    {
        const word& phase = thermoData_.phase(species[s]);
        data<< " =========================================================="
            << "=========================================================\n"
            << " Thermo analyses for " << species[s] << "\n"
            << " =========================================================="
            << "=========================================================\n"
            << " Phase: " << phase << "\n"
            << " Formula: " << formula[s] << "\n"
            << " Composition: \n   |\n";

        const map<word, scalar>& atoms =
            thermoData_.atomsAndFactors(species[s]);

        map<word, scalar> t = atoms;

        forMap(t, a)
        {
            data<< "   |--> " << std::setw(4) << a->first
                << " " << a->second << "\n";
        }

        data<< "\n" << std::setw(40) << std::left
            <<" Molecular weight:   "
            << std::setw(10) << std::right
            << MW(species[s]) << " [g/mol]\n"
            << std::setw(40) << std::left
            << " Formation enthalpy (298K):   "
            << std::setw(10) << std::right
            << Hf(species[s]) * M << " [J/mol]\n"
            << std::setw(40) << std::left
            <<" Frormation free Gibbs energy (298K):   "
            << std::setw(10) << std::right
            << Gf(species[s]) * M << " [J/mol]\n";

        data<< std::right << "\n";

        data<< "\n\n ----------------------------------------------------------"
            << "---------------------------------------------------------\n"
            << "      T    |         cp              H               S       "
            << "        G              dH             dG           | \n"
            << "     [K]   |      [J/molK]        [J/mol]         [J/molK]   "
            << "     [J/mol]         [J/molK]       [J/mol]        | \n"
            << " ----------------------------------------------------------"
            << "---------------------------------------------------------\n";

        for(int i=300; i<=3000; i+=100)
        {
            data<< "  " << std::setw(6) << i << "   |" 
                << "  " << std::setw(13) << cp(species[s], i) * M
                << "  " << std::setw(14) << H(species[s], i)  * M
                << "  " << std::setw(14) << S(species[s], i)  * M
                << "  " << std::setw(14) << G(species[s], i)  * M
                << "  " << std::setw(14) << dH(species[s], i) * M
                << "  " << std::setw(14) << dG(species[s], i) * M
                << "     |\n";
        }

        data<< " ----------------------------------------------------------"
            << "---------------------------------------------------------\n\n"
            << "\n\n";
    }
}


// ************************************************************************* //
