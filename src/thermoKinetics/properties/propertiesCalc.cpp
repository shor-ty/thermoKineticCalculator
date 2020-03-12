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

AFC::PropertiesCalc::PropertiesCalc
(
    const string fileName,
    const Thermo& thermo,
    const Chemistry& chemistry
)
:
    PropertiesData(fileName),

    thermo_(thermo),

    chemistry_(chemistry)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::PropertiesCalc::~PropertiesCalc()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
AFC::scalar AFC::PropertiesCalc::Omin
(
    const map<word, scalar>& fuelZj,
    const Thermo& thermo
)
{
    const scalar OminC = thermo.MW("O2")/thermo.MW("C");
    const scalar OminH = thermo.MW("O2")/(2*thermo.MW("H2"));

    //- Minimum oxigen
    scalar Omin{0};

    loopMap(species, value, fuelZj)
    {
        if (species == "C")
        {
            Omin += value*OminC;
        }
        else if (species == "H")
        {
            Omin += value*OminH;
        }
        else if (species == "O")
        {
            Omin -= value;
        }
    } 

    return Omin;
}


AFC::scalar AFC::PropertiesCalc::Zst
(
    const map<word, scalar>& fuelZj,
    const scalar Yf,
    const scalar YO2oxidizer,
    const Thermo& thermo
)
{
    const scalar omin = PropertiesCalc::Omin(fuelZj, thermo);

    Info<< "omin = "<< omin << endl;
    Info<< "Yf   = "<< Yf << endl;
    Info<< "YO2  = "<< YO2oxidizer << endl;
    return pow(1+(Yf*omin)/YO2oxidizer, -1);
}


AFC::scalar AFC::PropertiesCalc::adiabateFlameTemperature
(
    const PropertiesCalc& prop,
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

            //- First 'p' is an element in any case
            if (i == 1)
            {
                element = p;
                speciesComposition[s][element] += 0;
            } 
            //- Last 'p' has a special treatment
            else if (i == nCharacters)
            {
                //- If it is a digit add the atom character and add everything
                //  to the last element
                if (std::isdigit(p))
                {
                    //- Add character to string
                    nAtoms += p;

                    //- If last character, add the atoms to the last element
                    if (i == nCharacters)
                    {
                        speciesComposition[s][element] += std::stoi(nAtoms);
                    }
                }
                //- If it is not a digit, it is a single element at the end
                else
                {
                    //- Add one atom to previous element if the previous
                    //  character was an element too
                    if (nAtoms.empty())
                    {
                        speciesComposition[s][element] += 1;
                    }

                    element = p;
                    speciesComposition[s][element] += 1;
                }
            }
            //- Second and next p's are atoms and elements
            else
            {
                if (std::isdigit(p))
                {
                    //- Add character to string
                    nAtoms += p;

                    //- If last character, add the atoms to the last element
                    if (i == nCharacters)
                    {
                        speciesComposition[s][element] += std::stoi(nAtoms);
                    }
                }
                else
                {
                    //- If atoms is empty, the previous loop was also an
                    //  element so add one atom
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


AFC::map<AFC::word, AFC::scalar> AFC::PropertiesCalc::elementMassFraction
(
    const map<word, map<word, unsigned int>>& speciesComposition,
    const map<word, scalar>& Y,
    const Thermo& thermo
)
{
    wordList elements;

    //- Get all elements
    loopMap(species, composition, speciesComposition)
    {
        loopMap(element, nAtoms, composition)
        {
            //- Check if already inserted
            if (!elements.empty())
            {
                bool found{false};

                forAll(elements, e)
                {
                    if (e == element)
                    {
                        found = true;
                        break;
                    }    
                }

                //- No duplication
                if (!found)
                {
                    elements.push_back(element);
                }
            }
            //- First entry
            else
            {
                elements.push_back(element);
            }
        }
    }

    //- Calculate the mass fraction of the elements
    map<word, scalar> Zj;
    unsigned int i{0};

    forAll(elements, element)
    {
        Zj[element] = 0;

        loopMap(species, composition, speciesComposition)
        {
            const scalar speciesMW = thermo.MW(species);
            const scalar elementMW = thermo.MW(element);

            //- Copy as otherwise no 'find' fucntion can be used
            //  TODO should work in any case
            map<word, unsigned int> tmp = composition;
            map<word, unsigned int>::iterator it = tmp.find(element);
            unsigned int nAtoms{0};

            //- Element not found, set zero
            if (it != tmp.end())
            {
                nAtoms = it->second;
            }

            Zj[element] += nAtoms*elementMW*Y.at(species)/speciesMW;
        }
    }

    return Zj;
}


AFC::map<AFC::word, AFC::scalar> AFC::PropertiesCalc::elementMolFraction
(
    const map<word, scalar>& Zj,
    const Thermo& thermo
)
{
    scalar M{0};

    //- Calculate the Molecular weight
    loopMap(element, value, Zj)
    {
        const scalar elementMW = thermo.MW(element);

        M += value/elementMW;
    }

    M = 1/M;

    //- Calculate the mole fraction of the elements
    map<word, scalar> Wj;
    
    loopMap(element, value, Zj)
    {
        Wj[element] = value * M / thermo.MW(element);
    }

    return Wj;
}
*/

// ************************************************************************* //
