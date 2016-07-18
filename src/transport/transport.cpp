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

#include "transport.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Transport::Transport
(
    const string& fileName,
    const Thermo& thermo
)
    : thermo_(thermo)
{
    TransportReader transReader(fileName);
    
    transReader.read(transData_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Transport::~Transport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AFC::Transport::insertChemistrySpecies
(
    const wordList& speciesCH 
)
{
    transData_.insertChemistrySpecies(speciesCH);
}


void AFC::Transport::prepareFitting
(
    const Thermo& thermo 
)
{
    Info<< " c-o Prepare values for fitting procedure\n" << endl;

    //- Calculate viscosity using gas kinetics
    /*transCalc_.viscosity(thermo, transData_);

    //- Calculate binary diffusivity using gas kinetics
    transCalc_.binaryDiffusivity(thermo, transData_);

    //- Calculate thermal conductivity using gas kinetics
    transCalc_.thermalConductivity(thermo, transData_);
    */
}


AFC::scalar AFC::Transport::viscosity
(
    const word& species,
    const scalar& T,
    const word& method
) const
{
    return 
    (
        transCalc_.viscosity
        (
            species,
            T,
            thermo_,
            transData_,
            method
        )
    );
}


AFC::scalar AFC::Transport::thermalConductivity
(
    const word& species,
    const scalar& T,
    const word& method
) const
{
    return
    (
        transCalc_.thermalConductivity
        (
            species,
            T,
            thermo_,
            transData_,
            method
        )
    );
}


AFC::scalar AFC::Transport::binaryDiffusivity
(
    const word& species1,
    const word& species2,
    const scalar& T,
    const word& method
) const
{
    return 
    (
        transCalc_.binaryDiffusivity
        (
            species1,
            species2,
            T,
            thermo_,
            transData_,
            method
        )
    );
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

AFC::wordList AFC::Transport::species() const
{
    return transData_.species();
}


AFC::wordList AFC::Transport::chemistrySpecies() const
{
    return transData_.chemistrySpecies();
}


AFC::wordList AFC::Transport::chemicalFormula() const
{
    return transData_.chemicalFormula();
}


AFC::word AFC::Transport::chemicalFormula
(
    const word& species 
) const
{
    return transData_.chemicalFormula(species);
}


int AFC::Transport::geometricalConfig
(
    const word& species
) const
{
    return transData_.geometricalConfig(species);
}


AFC::scalar AFC::Transport::LJCD
(
    const word& species
) const
{
    return transData_.LJCD(species);
}


AFC::scalar AFC::Transport::LJP
(
    const word& species
) const
{
    return transData_.LJP(species);
}


AFC::scalar AFC::Transport::muk
(
    const word& species
) const
{
    return transData_.muk(species);
}


AFC::scalar AFC::Transport::alpha
(
    const word& species
) const
{
    return transData_.alpha(species);
}


AFC::scalar AFC::Transport::ZRot298
(
    const word& species
) const
{
    return transData_.ZRot298(species);
}


// * * * * * * * * * * * * * * * Summary Functions * * * * * * * * * * * * * //

void AFC::Transport::summary
(
    ostream& data 
) const
{
    data<< Header() << "\n";

    data<< " c-o Transport summary:\n"
        << "=======================\n\n"
        << "======================================================"
        << "===================================================\n"
        << "       Species      |   Geo. Config.    eps/kb    "
        << "    sigma          mu        alpha      zRot293    |\n"
        << "======================================================"
        << "===================================================\n";

    //  TODO switch to plot all or onyl the needed
    //- Just use the species that are used in chemistry
    const wordList& species_ = chemistrySpecies();

    data<<std::fixed;
    data.precision(4);

    forAll(species_, s)
    {
        data<< "  " << std::left << std::setw(15) << s << "   |"
            << std::right
            << "  " << std::setw(8) << geometricalConfig(s)
            << "  " << std::setw(13) << LJCD(s)
            << "  " << std::setw(12) << LJP(s)
            << "  " << std::setw(10) << muk(s)
            << "  " << std::setw(10) << alpha(s)
            << "  " << std::setw(10) << ZRot298(s)
            << "     |\n";
    }

    data<< "======================================================"
        << "===================================================\n";

    data<< "\n\n";

    forAll(species_, s)
    {
        data<< "============================================================="
            << "==\n"
            << " c-o Transport analyses for " << s << "\n"
            << "============================================================="
            << "==\n";

        data<< "\n\n"
            << "-------------------------------------------------------------"
            << "--\n"
            << "      T    |         mu           lambda         Dij       |\n"
            << "     [K]   |       [Pa s]        [W/m/K]        [J/molK]   |\n"
            << "-------------------------------------------------------------"
            << "--\n";

        data<<std::scientific;

        for(int i=300; i<=3000; i+=100)
        {
            data<< "  " << std::setw(6) << i << "   |" 
                << "  " << std::setw(13) << viscosity(s, i) 
                << "  " << std::setw(13) << thermalConductivity(s, i)
                << "  " << std::setw(13) << binaryDiffusivity("HE", "AR", i)
                << "     |\n";
        }

        data<< "-------------------------------------------------------------\n\n"
            << "\n\n";
    }
}


// ************************************************************************* //
