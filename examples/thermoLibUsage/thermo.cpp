/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2015 Holzmann-cfd
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Creator.

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

Description
    Example how to use the AFC modules for using the thermodynamic library
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
    AFC::Thermo thermo("NASA");

    //- Get the species for which we do have the NASA polynomials
    const AFC::wordList& thermoSpecies = thermo.species();

    AFC::Info<< "The NASA polynomials file does provide " 
        << thermoSpecies.size() << " species namely:\n\n";

    //- c++11 initialization
    int i{0};

    //- AFC loop usage
    for(auto species : thermoSpecies)
    {
        AFC::Info<< std::setw(18) << species;

        if (i == 5)
        {
            AFC::Info<< "\n";
            i = -1;
        }
        ++i;
    }

    AFC::Info<< "\n\nThe molecular weight of H2O is: " << thermo.MW("H2O")
        << " [g/mol]" << AFC::endl;

    AFC::Info<< "The molecular weight of H2H4 is: " << thermo.MW("N2H4")
        << " [g/mol]" << AFC::endl;

    AFC::Info<< "The molecular weight of C10H7CHO is: " << thermo.MW("C10H7CHO")
        << " [g/mol]" << AFC::endl;

    AFC::Info<< "The specie C4H9OOH consists of the following elements:"
        << AFC::endl;

    const AFC::wordList& elements = thermo.elementsInSpecies("C4H9OOH");
    const AFC::scalarList& atoms = thermo.elementAtoms("C4H9OOH");

    for (unsigned int i = 0; i < elements.size(); ++i)
    {
        AFC::Info<< "--> " << atoms[i] << " x " << elements[i] << "\n";
    }

    const AFC::map<AFC::word, AFC::scalar> elementAndAtoms =
        thermo.elementAtomsMap("C10H7CHO");


    AFC::Info<< "The specie C10H7CHO consists of the following elements:"
        << AFC::endl;

    //- Easy map loop
    for(auto [element, atoms] : elementAndAtoms)
    {
        AFC::Info<< element << " x " << atoms << "\n";
    }

    AFC::Info<< "The reaction enthalpy of H2O(504K) = "
        << thermo.h("H2O", 504) << " [J/mol]" << AFC::endl;

    AFC::Info<< "The heat capacity of NH3(292K) = "
        << thermo.cp("NH3", 292) << " [J/mol/K]" << AFC::endl;

    AFC::Info<< "The entropy of H2(1000K) = "
        << thermo.s("N2", 1000) << " [J/mol]" << AFC::endl;



    //- Bring output to terminal
    AFC::Info<< AFC::endl;

    return 0;
}


// ************************************************************************* //
