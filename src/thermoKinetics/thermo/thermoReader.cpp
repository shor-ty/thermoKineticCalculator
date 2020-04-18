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

#include "thermoReader.hpp" 
#include "constants.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::ThermoReader::ThermoReader(const string file, ThermoData& data)
:
    file_(file),
    data_(data)
{
    if (debug_)
    {
        Info<< "ThermoReader Constructor\n" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::ThermoReader::~ThermoReader()
{
    if (debug_)
    {
        Info<< "ThermoReader Destructor\n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AFC::ThermoReader::read()
{
    Info<< " c-o Reading the thermodynamic data (NASA Polynomials)\n"
        << "     >> " << file_ << "\n" << endl;

    const auto fileContent = readFile(file_);

    int lineNoKeyword{-1};
    unsigned int lineNoEnd{0};

    findKeyword(lineNoKeyword, lineNoEnd, fileContent);

    //- Reading THERMO block
    for (unsigned int line = lineNoKeyword+2; line < lineNoEnd; line++)
    {
        string ttmp = fileContent[line];

        //- Remove any comment
        removeComment(ttmp);

        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if (!tmp.empty())
        {
            if (fileContent[line][79] == '1')
            {
                setNASAPolynomialCoeffsNo1(fileContent[line], line);
                setNASAPolynomialCoeffsNo2(fileContent[++line], line);
                setNASAPolynomialCoeffsNo3(fileContent[++line], line);
                setNASAPolynomialCoeffsNo4(fileContent[++line], line);
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

    forEach(fileContent, line)
    {
        //- Split string; delimiter ' ' 
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- Search line no.
        forEach(searchPattern, i)
        {
            if (!tmp.empty() && tmp[0] == searchPattern[i])
            {
                start = line; 
                break;
            }
        }

        //- Search end after keyword found
        if (!tmp.empty() && tmp[0] == "END" && start != -1)
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
    const word element,
    const scalar multiplicator,
    const word species 
)
{
    if (!(AFC::Constants::AW.find(element) != AFC::Constants::AW.end()))
    {
        ErrorMsg
        (
            "    Element '" + element + "' not found in the data_base of the\n"
            "    thermoReader:: class. Please implement the element and its\n"
            "    weight into the thermoReader.hpp file. \n",
            __FILE__,
            __LINE__
        );
    }

    //- Multiplicator
    if (multiplicator <= 0)
    {
        ErrorMsg
        (
            "    Multiplicator of '" + element + " in species "
            + species + " is less or equal to zero.",
            __FILE__,
            __LINE__
        );
    }

    //- Element (atom) weight is in [g/mol]
    return AFC::Constants::AW.at(element)/scalar(1000) * multiplicator;
}


AFC::word AFC::ThermoReader::constructFormula(const string composition)
{
    //- Remove whitespaces
    wordList tmp = splitStrAtWS(composition);

    word formula;

    //- Form formula
    forEach(tmp, i)
    {
        formula += tmp[i];
    }

    return formula;
}


// * * * * * * * * * * * * Data manipulation functions * * * * * * * * * * * //

void AFC::ThermoReader::setNASAPolynomialCoeffsNo1
(
    const string lineContent,
    const unsigned int& line
)
{
    // Species name [1-18]
    wordList species = splitStrAtWS(lineContent.substr(0,18));
    {
        if (!species[0].empty())
        {
            data_.setSpecies(species[0]);
        }
        else
        {
            ErrorMsg
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

        //- Get all information about the elements and factors of the species
        const word elementicComposition = lineContent.substr(24,20);

        //- Build the species (formula)
        const word formula = constructFormula(elementicComposition);

        //- Insert chemical formula
        data_.setChemicalFormula(formula);

        //- Split formula into elements and factors and store them
        elementsAndFactors(formula);

        //- Calculate molecular weight and store them
        calcMolecularWeight(species[0]);
    }

    //- Phase [45]
    data_.setPhase(lineContent.substr(44,1));

    //- Low temperature [46-55]
    data_.setLT(stod(lineContent.substr(45,10)));

    //- High temperature [56-65]
    data_.setHT(stod(lineContent.substr(55,10)));

    //- Common temperature [66-73]
    data_.setCT(stod(lineContent.substr(65,8)));

    //- Addition elementic symbolic and formula [74-78] 
    //  Not implemented yet (should be done before)

    //- Integer 1 [80]

    //- Addition elementic symbolic and formula [81-100]
    //  Not implemented yet (should be done before)
}


void AFC::ThermoReader::setNASAPolynomialCoeffsNo2
(
    const string lineContent,
    const unsigned int& line
)
{
    //- First check integer '2' [80]
    if (lineContent[79] != '2')
    {
        ErrorMsg
        (
            "    Thermodynamic data_base is destroyed in line "
            + std::to_string(line) + ". No. '2' not found at pos 80.",
            __FILE__,
            __LINE__
        );
    }
    
    //- Coefficient a1 [1-15]
    data_.setNASACoeffsHT(stod(lineContent.substr(0,15)));

    //- Coefficient a2 [16-30]
    data_.setNASACoeffsHT(stod(lineContent.substr(15,15)));

    //- Coefficient a3 [31-45]
    data_.setNASACoeffsHT(stod(lineContent.substr(30,15)));

    //- Coefficient a4 [46-60]
    data_.setNASACoeffsHT(stod(lineContent.substr(45,15)));

    //- Coefficient a5 [61-75]
    data_.setNASACoeffsHT(stod(lineContent.substr(60,15)));
}


void AFC::ThermoReader::setNASAPolynomialCoeffsNo3
(
    const string lineContent,
    const unsigned int& line
)
{
    //- First check integer '3' [80]
    if (lineContent[79] != '3')
    {
        ErrorMsg
        (
            "    Thermodynamic data_base is destroyed in line "
            + std::to_string(line) + ". No. '3' not found at pos 80.",
            __FILE__,
            __LINE__
        );
    }
    
    //- Coefficient a6 [1-15]
    data_.setNASACoeffsHT(stod(lineContent.substr(0,15)));

    //- Coefficient a7 [16-30]
    data_.setNASACoeffsHT(stod(lineContent.substr(15,15)));

    //- Coefficient b1 [31-45]
    data_.setNASACoeffsLT(stod(lineContent.substr(30,15)));

    //- Coefficient b2 [46-60]
    data_.setNASACoeffsLT(stod(lineContent.substr(45,15)));

    //- Coefficient b3 [61-75]
    data_.setNASACoeffsLT(stod(lineContent.substr(60,15)));
}


void AFC::ThermoReader::setNASAPolynomialCoeffsNo4
(
    const string lineContent,
    const unsigned int& line
)
{
    //- First check integer '4' [80]
    if (lineContent[79] != '4')
    {
        ErrorMsg
        (
            "    Thermodynamic data_base is destroyed in line "
            + std::to_string(line) + ". No. '4' not found at pos 80.",
            __FILE__,
            __LINE__
        );
    }
    
    //- Coefficient b4 [1-15]
    data_.setNASACoeffsLT(stod(lineContent.substr(0,15)));

    //- Coefficient b5 [16-30]
    data_.setNASACoeffsLT(stod(lineContent.substr(15,15)));

    //- Coefficient b6 [31-45]
    data_.setNASACoeffsLT(stod(lineContent.substr(30,15)));

    //- Coefficient b7 [46-60]
    data_.setNASACoeffsLT(stod(lineContent.substr(45,15)));
}


void AFC::ThermoReader::calcMolecularWeight(const word species)
{
    //- Elements of species
    const wordList& elements = data_.elementsInSpecies(species);

    //- Factors of elements 
    const map<word, scalar>& factors = data_.elementAtomsMap(species);

    //- Temp molecular weight
    scalar tmp = 0;

    forAll(elements, e)
    {
        tmp += calcWeight(e, factors.at(e), species);
    }

    //- Molecular weight [kg/mol]
    data_.setMolecularWeight(tmp);
}


void AFC::ThermoReader::elementsAndFactors(const word formula)
{
    bool foundLetter{false};
    bool foundNumber{false};

    bool lastNumber{false};

    word element;

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
            //- Insert element and its factor
            data_.setElementAndAtoms(element, stoi(multiplicator));

            //- Reset
            element.clear();
            multiplicator.clear();
            lastNumber = false;
        }

        //- Form element
        if (foundLetter)
        {
            element += formula[pos];
        }
        else if (foundNumber)
        {
            lastNumber = true;

            multiplicator += formula[pos];
        }

        //- If pos is at last position
        if (pos == formula.size()-1)
        {
            //- Insert element and its factor
            data_.setElementAndAtoms(element, stoi(multiplicator));

            //- Everything is done
            break;
        }
    }

    //- All elements and factors stored, now we can update the map
    data_.updateElementsAndFactors();
}


// ************************************************************************* //
