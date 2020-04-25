/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     T hermo       | Open Source Thermo-Kinetic Library
  c-o-o-c     K iknetic     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Creator.

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

Description
    Example how to use the TKC modules for using the thermodynamic library
    Here we calculate the reaction enthalpy of different species out of the
    box as well as the free Gibbs energy


\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "definitions.hpp"
#include "thermo.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main ()
{

    //- Create the thermo object
    TKC::Thermo thermo("NASA");

    //- Get the species for which we do have the NASA polynomials
    const TKC::wordList& thermoSpecies = thermo.species();

    TKC::Info<< "The NASA polynomials file does provide " 
        << thermoSpecies.size() << " species namely:\n\n";

    //- c++11 initialization
    int i{0};

    //- TKC loop usage
    for(auto species : thermoSpecies)
    {
        TKC::Info<< std::setw(18) << species;

        if (i == 5)
        {
            TKC::Info<< "\n";
            i = -1;
        }
        ++i;
    }

    TKC::Info<< "\n\nThe molecular weight of H2O is: " << thermo.MW("H2O")
        << " [g/mol]" << TKC::endl;

    TKC::Info<< "The molecular weight of H2H4 is: " << thermo.MW("N2H4")
        << " [g/mol]" << TKC::endl;

    TKC::Info<< "The molecular weight of C10H7CHO is: " << thermo.MW("C10H7CHO")
        << " [g/mol]" << TKC::endl;

    TKC::Info<< "The specie C4H9OOH consists of the following elements:"
        << TKC::endl;

    const TKC::wordList& elements = thermo.elementsInSpecies("C4H9OOH");
    const TKC::scalarList& atoms = thermo.elementAtoms("C4H9OOH");

    for (unsigned int i = 0; i < elements.size(); ++i)
    {
        TKC::Info<< "--> " << atoms[i] << " x " << elements[i] << "\n";
    }

    const TKC::map<TKC::word, TKC::scalar> elementAndAtoms =
        thermo.elementAtomsMap("C10H7CHO");


    TKC::Info<< "The specie C10H7CHO consists of the following elements:"
        << TKC::endl;

    //- Easy map loop
    for(auto [element, atoms] : elementAndAtoms)
    {
        TKC::Info<< element << " x " << atoms << "\n";
    }

    TKC::Info<< "The reaction enthalpy of H2O(504K) = "
        << thermo.h("H2O", 504) << " [J/mol]" << TKC::endl;

    TKC::Info<< "The heat capacity of NH3(292K) = "
        << thermo.cp("NH3", 292) << " [J/mol/K]" << TKC::endl;

    TKC::Info<< "The entropy of H2(1000K) = "
        << thermo.s("N2", 1000) << " [J/mol]" << TKC::endl;



    //- Bring output to terminal
    TKC::Info<< TKC::endl;

    return 0;
}


// ************************************************************************* //
