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

#include "thermoReader.hpp" 
#include "constants.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::ThermoReader::ThermoReader
(
    const string& file
)
:
    file_(file)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::ThermoReader::~ThermoReader()
{
    Info<< "Destructor ThermoReader\n" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AFC::ThermoReader::read
(
    ThermoData& data
)
{
    Info<< " c-o Reading thermodynamic data\n" << endl;

    //- TODO if thermo true --> use chemistry file
    
    const auto fileContent = readFile(file_);

    int lineNoKeyword{-1};
    unsigned int lineNoEnd{0};

    findKeyword
    (
        lineNoKeyword,
        lineNoEnd,
        fileContent
    );

    //- Reading THERMO block
    for (unsigned int line = lineNoKeyword+2; line < lineNoEnd; line++)
    {
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            if (fileContent[line][79] == '1')
            {
                NASAPolynomialNo1(fileContent[line], line, data);
                NASAPolynomialNo2(fileContent[++line], line, data);
                NASAPolynomialNo3(fileContent[++line], line, data);
                NASAPolynomialNo4(fileContent[++line], line, data);
            }
        }
    }
}


// * * * * * * * * * * * * * * Helper functions  * * * * * * * * * * * * * * //

void AFC::ThermoReader::findKeyword
(
    int& start,
    unsigned int& end,
    const stringList& fileContent
)
{
    wordList searchPattern = THERMO;

    forAll(fileContent, line)
    {
        //- Split string; delimiter ' ' 
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- Search line no.
        forAll(searchPattern, i)
        {
            if (tmp[0] == searchPattern[i])
            {
                start = line; 
                break;
            }
        }

        //- Search end after keyword found
        if (tmp[0] == "END" && start != -1)
        {
            end = line;
        }

        //- Both found exit loop
        if
        (
            start != 0
         && end != 0
        )
        {
            break;
        }
    }
}


AFC::scalar AFC::ThermoReader::calcWeight
(
    const word& atom,
    const scalar& multiplicator,
    const word& species 
)
{
    if (!(AFC::Constants::AW.find(atom) != AFC::Constants::AW.end()))
    {
        FatalError
        (
            "    Atom '" + atom + "' not found in the database of the\n"
            "    thermoReader:: class. Please implement the atom and its\n"
            "    weight into the thermoReader.hpp file. \n",
            __FILE__,
            __LINE__
        );
    }

    //- Multiplicator
    if (multiplicator <= 0)
    {
        FatalError
        (
            "    Multiplicator is less than 0 of '" + atom + " in species "
            + species + ".",
            __FILE__,
            __LINE__
        );
    }

    return AFC::Constants::AW.at(atom) * multiplicator;
}


AFC::word AFC::ThermoReader::constructFormula
(
    const string& composition 
)
{
    //- Remove whitespaces
    wordList tmp = splitStrAtWS(composition);

    word formula;

    //- Form formula
    forAll(tmp, i)
    {
        formula += tmp[i];
    }

    return formula;
}


// * * * * * * * * * * * * Data manipulation functions * * * * * * * * * * * //

void AFC::ThermoReader::NASAPolynomialNo1
(
    const string& lineContent,
    const unsigned int& line,
    ThermoData& data
)
{
    // Species name [1-18]
    wordList species = splitStrAtWS(lineContent.substr(0,18));
    {
        if (!species[0].empty())
        {
            data.insertSpecies(species[0]);
        }
        else
        {
            FatalError
            (
                "    No species entry found in file " + file_ + " line no. "
                + std::to_string(line),
                __FILE__,
                __LINE__
            );
        }
    }

    // Species formula [22-44] || more complex stuff
    {

        //- Get all information about the atoms and factors of the species
        const word atomicComposition = lineContent.substr(24,20);

        //- Build the species (formula)
        const word formula = constructFormula(atomicComposition);

        //- Insert chemical formula
        data.insertChemicalFormula(formula);

        //- Split formula into atoms and factors and store them
        atomsAndFactors(formula, data);

        //- Calculate molecular weight and store them
        calcMolecularWeight(species[0], data);
    }

    //- Phase [45]
    data.insertPhase(lineContent.substr(44,1));

    //- Low temperature [46-55]
    data.insertLT(stod(lineContent.substr(45,10)));

    //- High temperature [56-65]
    data.insertHT(stod(lineContent.substr(55,10)));

    //- Common temperature [66-73]
    data.insertCT(stod(lineContent.substr(65,8)));

    //- Addition atomic symbolic and formula [74-78] 
    //  Not implemented yet (should be done before)

    //- Integer 1 [80]

    //- Addition atomic symbolic and formula [81-100]
    //  Not implemented yet (should be done before)
}


void AFC::ThermoReader::NASAPolynomialNo2
(
    const string& lineContent,
    const unsigned int& line,
    ThermoData& data
)
{
    //- First check integer '2' [80]
    if (lineContent[79] != '2')
    {
        FatalError
        (
            "    Thermodynamic database is destroyed in line "
            + std::to_string(line) + ". No. '2' not found at pos 80.",
            __FILE__,
            __LINE__
        );
    }
    
    //- Coefficient a1 [1-15]
    data.insertNASACoeffsHT(stod(lineContent.substr(0,15)));

    //- Coefficient a2 [16-30]
    data.insertNASACoeffsHT(stod(lineContent.substr(15,15)));

    //- Coefficient a3 [31-45]
    data.insertNASACoeffsHT(stod(lineContent.substr(30,15)));

    //- Coefficient a4 [46-60]
    data.insertNASACoeffsHT(stod(lineContent.substr(45,15)));

    //- Coefficient a5 [61-75]
    data.insertNASACoeffsHT(stod(lineContent.substr(60,15)));
}


void AFC::ThermoReader::NASAPolynomialNo3
(
    const string& lineContent,
    const unsigned int& line,
    ThermoData& data
)
{
    //- First check integer '3' [80]
    if (lineContent[79] != '3')
    {
        FatalError
        (
            "    Thermodynamic database is destroyed in line "
            + std::to_string(line) + ". No. '3' not found at pos 80.",
            __FILE__,
            __LINE__
        );
    }
    
    //- Coefficient a6 [1-15]
    data.insertNASACoeffsHT(stod(lineContent.substr(0,15)));

    //- Coefficient a7 [16-30]
    data.insertNASACoeffsHT(stod(lineContent.substr(15,15)));

    //- Coefficient b1 [31-45]
    data.insertNASACoeffsLT(stod(lineContent.substr(30,15)));

    //- Coefficient b2 [46-60]
    data.insertNASACoeffsLT(stod(lineContent.substr(45,15)));

    //- Coefficient b3 [61-75]
    data.insertNASACoeffsLT(stod(lineContent.substr(60,15)));
}


void AFC::ThermoReader::NASAPolynomialNo4
(
    const string& lineContent,
    const unsigned int& line,
    ThermoData& data
)
{
    //- First check integer '4' [80]
    if (lineContent[79] != '4')
    {
        FatalError
        (
            "    Thermodynamic database is destroyed in line "
            + std::to_string(line) + ". No. '4' not found at pos 80.",
            __FILE__,
            __LINE__
        );
    }
    
    //- Coefficient b4 [1-15]
    data.insertNASACoeffsLT(stod(lineContent.substr(0,15)));

    //- Coefficient b5 [16-30]
    data.insertNASACoeffsLT(stod(lineContent.substr(15,15)));

    //- Coefficient b6 [31-45]
    data.insertNASACoeffsLT(stod(lineContent.substr(30,15)));

    //- Coefficient b7 [46-60]
    data.insertNASACoeffsLT(stod(lineContent.substr(45,15)));
}


void AFC::ThermoReader::calcMolecularWeight
(
    const word& species,
    ThermoData& data
)
{
    //- Atoms of species
    const wordList& atoms = data.speciesAtoms(species);

    //- Factors of atoms
    const scalarList& factors = data.atomFactors(species);

    //- Temp molecular weight
    scalar tmp = 0;

    forAll(atoms, a)
    {
        tmp += calcWeight(atoms[a], factors[a], species);
    }

    //- Molecular weight [g/mol]
    data.insertMolecularWeight(tmp);
}


void AFC::ThermoReader::atomsAndFactors
(
    const word& formula,
    ThermoData& data
)
{
    bool foundLetter{false};
    bool foundNumber{false};

    bool lastNumber{false};

    word atom;

    //- Type string due to adding letters
    string multiplicator{""};
    scalar tmpMW{0};

    //- Loop through all single letter of the formula (species)
    for (unsigned int pos = 0; pos < formula.size(); pos++)
    {
        foundLetter = false;
        foundNumber = false;

        //- Compare single letter with ASCII table (LETTERS)
        for (unsigned int j = 65; j <= 90; j++)
        {
            char c = static_cast<char>(j);

            if
            (
                c == formula[pos]
            )
            {
                foundLetter = true; foundNumber = false;
                break;
            }
            else
            {
                foundLetter = false;
                foundNumber = true;
            }
        }

        //- If letter found and last letter was a number (new element)
        if (foundLetter && lastNumber)
        {
            //- Insert atom and its factor
            data.insertAtomAndFactor(atom, stoi(multiplicator));

            //- Reset
            atom.clear();
            multiplicator.clear();
            lastNumber = false;
        }

        //- Form atom
        if (foundLetter)
        {
            atom += formula[pos];
        }
        else if (foundNumber)
        {
            lastNumber = true;

            multiplicator += formula[pos];
        }

        //- If pos is at last position
        if (pos == formula.size()-1)
        {
            //- Insert atom and its factor
            data.insertAtomAndFactor(atom, stoi(multiplicator));

            //- Everything is done
            break;
        }
    }
}


// ************************************************************************* //
