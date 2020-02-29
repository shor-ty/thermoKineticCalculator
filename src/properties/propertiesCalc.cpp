/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
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

#include "propertiesCalc.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::PropertiesCalc::PropertiesCalc()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::PropertiesCalc::~PropertiesCalc()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

AFC::scalar AFC::PropertiesCalc::Zst
(
    const map<word, scalar>& Yfuel,
    const map<word, scalar>& Yoxidizer
)
{
    Info<< Yfuel.at("H2") << endl;    
    Info<< Yfuel.at("H2") << endl;    


    std::terminate();
    return 0;
}

AFC::scalar AFC::PropertiesCalc::adiabateFlameTemperature
(
    const Properties& prop,
    const Thermo& thermo,
    const Chemistry& chem
)
{
//    const wordList& fuel = prop.fuel();
    return 0;
}


AFC::map<AFC::word, AFC::map<AFC::word, unsigned int>>
AFC::PropertiesCalc::elementDecomposition
(
    const wordList& species,
    const Thermo& thermo 
)
{
    map<word, map<word, unsigned>> speciesComposition;
    map<word, unsigned int> composition;

    forAll(species, s)
    {
        // Add species to map
        speciesComposition[s] = map<word, unsigned>();

        size_t i = 0;
        size_t nCharacters = s.size();
        string nAtoms = "";
        string element = "";
        
        //- Move through the species and find elements + amount of atoms
        //  Note: Structure is always 'Element + digits'
        forAll(s, p)
        {
            ++i; 

            //- First 'p' is first element
            if (i == 1)
            {
                element = p;
                speciesComposition[s][element] += 0;
            } 
            //- Second and next p's are atoms and elements
            else
            {
                if (std::isdigit(p))
                {
                    //- Increment
                    nAtoms += p;

                    //- If last character, add the atoms to the last element
                    if (i == nCharacters)
                    {
                        speciesComposition[s][element] += std::stoi(nAtoms);
                    }
                }
                else
                {
                    //- Futhermore, if nAtoms is zero, add one atom
                    if (nAtoms.empty())
                    {
                        speciesComposition[s][element] += 1;
                    }
                    else
                    {
                        speciesComposition[s][element] += std::stoi(nAtoms);
                        nAtoms = "";
                    }

                    //- Actual position is new element
                    element = p;
                }
            }
        }
    }

    return std::move(speciesComposition);
}


AFC::map<AFC::word, unsigned int> AFC::PropertiesCalc::elementMassFraction
(
    const map<word, unsigned int>& composition,
    const Thermo& thermo
)
{
    loopMap(species, nAtoms, composition)
    {

    }
}


// ************************************************************************* //
