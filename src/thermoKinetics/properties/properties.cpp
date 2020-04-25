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

#include "properties.hpp"
#include "constants.hpp"
#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::Properties::Properties
(
    const string fileName,
    const Thermo& thermo,
    const Chemistry& chemistry
)
:
    PropertiesCalc(fileName, thermo, chemistry)
{
    //- Add pressure to thermoData for better handling
    //thermo.p(this->p());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::Properties::~Properties()
{}


// * * * * * * * * * * * * * * Summary Function  * * * * * * * * * * * * * * //

void TKC::Properties::summary(ostream& data) const
{
    //- Header
    data<< Header() << "\n"; 

    data<< " c-o Properties summary (analysis of afcDict):\n"
        << " =============================================\n\n\n"
        << " =============================================================="
        << "===============\n"
        << " c-o Oxidizer information\n"
        << " =============================================================="
        << "===============\n"
        << "  |\n"
        << "  |--> Number of species:    " << speciesOxidizer_.size() << "\n"
        << "  |--> Oxidizer set to be:   " << oxidizer_ << "\n"
        << "  | \n"
        << "  | "
        << std::setw(15) << "Species"
        << std::setw(20) << "Mass Fraction"
        << std::setw(20) << "Mole Fraction"
        << std::setw(20) << "Concentration\n"
        << "  | " 
        << std::setw(35) << " Y [-]"
        << std::setw(20) << " X [-]"
        << std::setw(20) << " C [g/mol]\n" 
        << "  |------------------------------------------------------------"
        << "---------------\n";

    scalar sumY{0};
    scalar sumX{0};
    forAll(speciesOxidizer_, species)
    {
        data<< "  | " << std::right << std::setw(15) << species; 
        data<< std::setw(20) << oxidizerY_.at(species);
        data<< std::setw(20) << oxidizerX_.at(species);
        //data<< std::setw(20) << oxidizerN_.at(species);
        data<< "\n";

        sumY += oxidizerY_.at(species);
        sumX += oxidizerX_.at(species);
    }

    data<< "  |------------------------------------------------------------"
        << "---------------\n"
        << "  | "<< std::setw(15) << "Sum"
        << std::setw(20) << sumY
        << std::setw(20) << sumX
        << std::setw(20) << "\n"
        << "  |------------------------------------------------------------"
        << "---------------\n";

    loopMap(element, value, oxidizerZj_)
    {
        data<< "  | " << std::right << std::setw(15) << element; 
        data<< std::setw(20) << value;
        data<< std::setw(20) << oxidizerWj_.at(element);
        //data<< std::setw(20) << oxidizerN_.at(species);
        data<< "\n";
    }

    data<< "  |------------------------------------------------------------"
        << "---------------\n";

    data<< "\n\n"
        << " =============================================================="
        << "===============\n"
        << " c-o Fuel information\n"
        << " =============================================================="
        << "===============\n"
        << "  |\n"
        << "  |--> Number of species:    " << speciesFuel_.size() << "\n"
        << "  |--> Fuel set to be:       " << fuel_ << "\n"
        << "  | \n"
        << "  |------------------------------------------------------------"
        << "---------------\n"
        << "  | "
        << std::setw(15) << "Species"
        << std::setw(20) << "Mass Fraction"
        << std::setw(20) << "Mole Fraction"
        << std::setw(20) << "Concentration\n"
        << "  | " 
        << std::setw(35) << " Y [-]"
        << std::setw(20) << " X [-]"
        << std::setw(20) << " C [g/mol]\n" 
        << "  |------------------------------------------------------------"
        << "---------------\n";

    sumY = 0;
    sumX = 0;
    forAll(speciesFuel_, species)
    {
        data<< "  | " << std::right << std::setw(15) << species; 
        data<< std::setw(20) << fuelY_.at(species);
        data<< std::setw(20) << fuelX_.at(species);
        //data<< std::setw(20) << oxidizerN_.at(species);
        data<< "\n";

        sumY += fuelY_.at(species);
        sumX += fuelX_.at(species);
    }

    data<< "  |------------------------------------------------------------"
        << "---------------\n"
        << "  | "<< std::setw(15) << "Sum"
        << std::setw(20) << sumY
        << std::setw(20) << sumX
        << std::setw(20) << "\n"
        << "  |------------------------------------------------------------"
        << "---------------\n";

    loopMap(element, value, fuelZj_)
    {
        data<< "  | " << std::right << std::setw(15) << element; 
        data<< std::setw(20) << value;
        data<< std::setw(20) << fuelWj_.at(element);
        //data<< std::setw(20) << oxidizerN_.at(species);
        data<< "\n";
    }

    data<< "  |------------------------------------------------------------"
        << "---------------\n";
}
*/

// ************************************************************************* //
