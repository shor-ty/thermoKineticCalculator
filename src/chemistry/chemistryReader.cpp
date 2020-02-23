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

#include "chemistryReader.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::ChemistryReader::ChemistryReader(const string file)
:
    file_(file)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::ChemistryReader::~ChemistryReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AFC::ChemistryReader::read(ChemistryData& data)
{
    Info<< " c-o Reading chemistry data\n" << endl;

    const auto fileContent = readFile(file_);

    readElementBlock(fileContent, data);

    readSpeciesBlock(fileContent, data);

    readThermoBlock(fileContent, data);

    readReactionBlock(fileContent, data);
}


void AFC::ChemistryReader::readElementBlock
(
    const stringList& fileContent,
    ChemistryData& data
)
{
    //- STEP 1: find line no. of wordList ELEMENTS and "END"
    int lineNoKeyword{-1};
    unsigned int lineNoEnd{0};

    findKeyword
    (
        lineNoKeyword,
        lineNoEnd,
        fileContent,
        "E"
    );

    //- STEP 2: check if wordList ELEMENT found
    if
    (
        lineNoKeyword == -1
    )
    {
        ErrorMsg
        (
            "Keyword in list 'ELEMENT' not found in chemistry file " + file_,
            __FILE__,
            __LINE__
        );
    }

    //- Reading ELEMENT block
    for (unsigned int line = lineNoKeyword+1; line < lineNoEnd; line++)
    {
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            //- If more elements in one line
            if (tmp.size() > 1)
            {
                forAll(tmp, element)
                {
                    data.elements(element);
                }
            }
            else
            {
                data.elements(tmp[0]);
            }
        }
    }
}


void AFC::ChemistryReader::readSpeciesBlock
(
   const stringList& fileContent,
   ChemistryData& data
)
{
    //- STEP 1: find line no. of wordList SPECIES and "END"
    int lineNoKeyword{-1};
    unsigned int lineNoEnd{0};

    findKeyword
    (
        lineNoKeyword,
        lineNoEnd,
        fileContent,
        "S"
    );

    //- STEP 2: check if wordList SPECIES found
    if
    (
        lineNoKeyword == -1
    )
    {
        ErrorMsg
        (
            "Keyword in list 'SPECIES' not found in chemistry file " + file_,
            __FILE__,
            __LINE__
        );
    }

    //- Reading SPECIES block
    for (unsigned int line = lineNoKeyword+1; line < lineNoEnd; line++)
    {
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            //- If more species in one line
            if (tmp.size() > 1)
            {
                forAll(tmp, species)
                {
                    data.species(species);
                }
            }
            else
            {
                data.species(tmp[0]);
            }
        }
    }
}


void AFC::ChemistryReader::readThermoBlock
(
    const stringList& fileContent,
    ChemistryData& data
)
{
    //- STEP 1: find line no. of wordList THERMO and "END"
    int lineNoKeyword{-1};
    unsigned int lineNoEnd{0};

    findKeyword
    (
        lineNoKeyword,
        lineNoEnd,
        fileContent,
        "T"
    );

    //- STEP 2: check if wordList THERMO found
    if (lineNoKeyword != -1)
    {
        data.setThermo();
    }
}

void AFC::ChemistryReader::readReactionBlock
(
    const stringList& fileContent,
    ChemistryData& data
)
{
    //- STEP 1: find line no. of wordList REACTIONS and "END"
    int lineNoKeyword{-1};
    unsigned int lineNoEnd{0};

    findKeyword
    (
        lineNoKeyword,
        lineNoEnd,
        fileContent,
        "R"
    );

    //- STEP 2: check if wordList REACTION found
    if (lineNoKeyword == -1)
    {
        ErrorMsg
        (
            "Keyword in list 'REACTIONS' not found in chemistry file " + file_,
            __FILE__,
            __LINE__
        );
    }

    //- Reading REACTION block
    for (unsigned int line = lineNoKeyword+1; line < lineNoEnd; line++)
    {
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and is not a comment, proceed
        if
        (
            !tmp.empty()
            && tmp[0][0] != '!'
        )
        {
            //- Check for duplicate entrys
            if (tmp[0] == "DUPLICATE")
            {
                //- Skip next two lines (third with for loop)
                line+= 2;

                //- Increment duplicate entry
                data.incrementDuplicated();
            }
            else
            {
                //- Check if another comment is somewhere in the line 
                std::size_t foundExMark = fileContent[line].find('!');

                if (foundExMark != std::string::npos)
                {
                    //- Search the entry of the comment
                    forEach(tmp, i)
                    {
                        if (tmp[i][0] == '!')
                        {
                            tmp.resize(i);
                        } 
                    };
                }

                //- Check if '=' is in string (means reaction)
                std::size_t found = fileContent[line].find('=');

                if (found != std::string::npos)
                {
                    data.incrementReac();

                    data.incrementMatrixesVectors();

                    analyzeReaction(fileContent[line], data);

                    found = fileContent[line].find("+M");

                    if (found != std::string::npos)
                    {
                        data.TBR(true);
                    }
                }
                else
                {
                    std::size_t foundLOW;
                    std::size_t foundSRI;
                    std::size_t foundTROE;

                    foundLOW = fileContent[line].find("LOW");
                    foundSRI = fileContent[line].find("SRI");
                    foundTROE = fileContent[line].find("TROE");

                    //- LOW parameters
                    if (foundLOW != std::string::npos)
                    {
                        data.LOW(true);
                        LOWCoeffs(fileContent[line], line, data);
                    }
                    //- TROE parameters
                    else if (foundTROE != std::string::npos)
                    {
                        data.TROE(true);
                        TROECoeffs(fileContent[line], line, data);
                    }
                    //- SRI parameters
                    else if (foundSRI != std::string::npos)
                    {
                        data.SRI(true);
                        SRICoeffs(fileContent[line], line, data);
                    }
                    else
                    {
                        data.ENHANCE(true);
                        enhanceFactors(fileContent[line], data);
                    }
                }
            }
        }
    }

    //- At the end we have to increase nReac because it is initialized as -1
    data.incrementReac();
}


// * * * * * * * * * * * * * * Helper functions  * * * * * * * * * * * * * * //

void AFC::ChemistryReader::findKeyword
(
    int& start,
    unsigned int& end,
    const stringList& fileContent,
    const string search
)
{
    //- Search pattern (wordList)
    wordList searchPattern;

    // Check which section is needed

        //- ELEMENTS
        if (search == "E")
        {
            searchPattern = ELEMENT;
        }
        //- SPECIES
        else if (search == "S")
        {
            searchPattern = SPECIES;
        }
        //- THERMODYNAMIC
        else if (search == "T")
        {
            searchPattern = THERMO;
        }
        //- CHEMISTRY
        else if (search == "R")
        {
            searchPattern = REACTION;
        }

    // Iterator due to new forAll() Range-for 
    unsigned int lineNo{0};

    forAll(fileContent, line)
    {
        //- Split string; delimiter ' ' 
        stringList tmp = splitStrAtWS(line);

        //- Search line no.
        forAll(searchPattern, i)
        {
            if (tmp[0] == i)
            {
                start = lineNo; 
                break;
            }
        }

        //- Search end after keyword found
        if (tmp[0] == "END" && start != -1)
        {
            end = lineNo;
        }

        //- Both found exit loop
        if
        (
            start != -1
         && end != 0
        )
        {
            break;
        }

        lineNo++;
    }
}


AFC::stringList AFC::ChemistryReader::extractData(const string str)
{
    //- STEP 1: find first '/' 
    string delimiter="/";

    std::size_t found = str.find(delimiter);

    //- STEP 2: remove all letters from 0 till pos of '/'
    string tmp = str.substr(found+1,str.size());

    //- STEP 3: find second '/'
    found = tmp.find(delimiter);

    //- STEP 4: remove all letters from pos of '/' till end
    string tmp2 = tmp.substr(0,found);

    return splitStrAtWS(tmp2);
}


// * * * * * * * * * * * * Data manipulation functions * * * * * * * * * * * //

void AFC::ChemistryReader::analyzeReaction
(
    const string reaction,
    ChemistryData& data
)
{
    //- STEP 1: manipulate string to get reaction
    stringList tmp = splitStrAtWS(reaction);

    const string& tmp2 = tmp[0];

    //- Re-arrange the string and remove arrhenius coeffs
    /*for (unsigned int i=0; i<tmp.size()-3; i++)
    {
        tmp2 += tmp[i];
    }
    */

    data.elementarReaction(tmp2);
    
    //- STEP 2: analyze reaction
    string delimiter1 = "=";
    string delimiter2 = "<";
    string delimiter3 = ">";

    //- Reactions definition
    //
    //  +  a)  A+A=B   (forward and backward)
    //  +  b)  A+A<=>B (forward and backward)
    //  +  c)  A+A=>B  (only forward reaction)
    //  +  d)  A+A<=B  (only backward reaction)
    //  d) is not implemented

    size_t found1 = tmp2.find(delimiter1);
    size_t found2 = tmp2.find(delimiter2);
    size_t found3 = tmp2.find(delimiter3);

    //- a)
    if
    (
        found1 != std::string::npos
     && found2 == std::string::npos
     && found3 == std::string::npos
    )
    {
        data.BR(true);
    }
    //- b)
    else if
    (
        found1 != std::string::npos
     && found2 != std::string::npos
     && found3 != std::string::npos
    )
    {
        data.BR(true);
    }
    //- c)
    else if
    (
        found1 != std::string::npos
     && found2 == std::string::npos
     && found3 != std::string::npos
    )
    {
        data.BR(false);
    }
    //- d)
    //  Not implemented, normally not used
    //  NOTE: 14.07.16, change product and educt and fine
    else if
    (
        found1 != std::string::npos
     && found2 != std::string::npos
     && found3 == std::string::npos
    )
    {
        ErrorMsg
        (
            "Only backward reaction is not implemented. If for some reason\n   "
            " you want to implement it, feel free to open a issue at\n    "
            "www.bitbucket.org/shorty or implement it yourself and "
            "make a pull request.",
            __FILE__,
            __LINE__
        );
    }

    // STEP 3: insert arrhenius coeffs
    data.arrheniusCoeffs
    (
        stod(tmp[1]),
        stod(tmp[2]),
        stod(tmp[3])
    );

    // STEP 4: get stochiometric values
    //  + product negativ
    //  + reactants positiv

    // a) split into educts and products
    stringList tmp3 = splitStrAtDelimiter(tmp2, '=');

    const string reac = tmp3[0];
    const string prod = tmp3[1];

    // b) Check if the first char is a number,
    //    if not nu = 1
    //    if its a number, check next letter for number
    {
        //- Reactant site
        analyzeReacSite(reac, "e", data);

        //- Product site
        analyzeReacSite(prod, "p", data);
    }

    //- All stochiometric coefficients for species stored. Hence we
    //  are able to calculate the global reaction order and exponent
    //  for calculation Kc
    data.updateGlobalReactionOrder();
    
}


void AFC::ChemistryReader::LOWCoeffs
(
    const string coeffStr,
    const unsigned int lineNo,
    ChemistryData& data
)
{
    //- STEP 1: get data inbetween '/'
    List<string> coeffs = extractData(coeffStr);

    //- STEP 2: check if 3 values are available
    if (coeffs.size() != 3)
    {
        ErrorMsg
        (
            "More or less than three arrhenius coeffs for LOW found.\n    "
            "Arrhenius coeffs: " + std::to_string(coeffs.size()) +
            " found in line " + std::to_string(lineNo) + ".\n    "
            "Problem occur in file " + file_,
            __FILE__,
            __LINE__
        );
    }

    //- STEP 3: update 
    unsigned int c{0};

    forAll(coeffs, i)
    {
        data.LOWCoeffs(stod(i), c);
        c++;
    }
}


void AFC::ChemistryReader::TROECoeffs
(
    const string coeffStr,
    const unsigned int lineNo,
    ChemistryData& data
)
{
    //- STEP 1: get data inbetween '/'
    stringList coeffs = extractData(coeffStr);

    //- STEP 2: check if more or less values
    if (coeffs.size() < 3 || coeffs.size() > 4)
    {
        ErrorMsg
        (
            "More than 4 or less than 3 TROE coeffs found.\n    "
            "TROE coeffs: " + std::to_string(coeffs.size()) +
            " found in line " + std::to_string(lineNo) + ".\n    "
            "Problem occur in file " + file_,
            __FILE__,
            __LINE__
        );
    }

    //- STEP 3: update
    unsigned int c{0};

    forAll(coeffs, i)
    {
         data.TROECoeffs(stod(i), c);
         c++;
    }
}


void AFC::ChemistryReader::SRICoeffs
(
    const string coeffStr,
    const unsigned int lineNo,
    ChemistryData& data
)
{
    //- STEP 1: get data inbetween '/'
    stringList coeffs = extractData(coeffStr);

    //- STEP 2: check if more or less values
    if (coeffs.size() < 4 || coeffs.size() > 5)
    {
        ErrorMsg
        (
            "More than 5 or less than 4 SRI coeffs found.\n    "
            "SRI coeffs: " + std::to_string(coeffs.size()) +
            " found in line " + std::to_string(lineNo) + ".\n    "
            "Problem occur in file " + file_,
            __FILE__,
            __LINE__
        );
    }

    //- STEP 3: update
    unsigned int c{0};

    forAll(coeffs, i)
    {
         data.SRICoeffs(stod(i), c);
    }
}


void AFC::ChemistryReader::enhanceFactors
(
    const string enhanceFactors,
    ChemistryData& data
)
{
    //- STEP 1: split string at whitespace and re-arrange
    stringList tmp = splitStrAtWS(enhanceFactors);

    string enhanced;

    forAll(tmp, i)
    {
        enhanced += i;
    }

    //- STEP 2: split enhanced string at delimiter '/' and save
    tmp = splitStrAtDelimiter(enhanced, '/');

    //- Count entries
    unsigned int c{0};
    word species;

    forAll(tmp, i)
    {
        //- Species name
        if (c == 0)
        {
            species = i; 
            c++;
        }
        else if (c == 1)
        {
            data.ENHANCEDCoeffs(species, stod(i));
            c = 0;
        }
    }
}


void AFC::ChemistryReader::analyzeReacSite
(
    const string reactionSite,
    const word site,
    ChemistryData& data
)
{
    //- First: remove (+M) if it is there
    string removed1 = removeAtEnd(reactionSite, "(+M)");

    //- Second: remove +M if it is there
    string tmp = removeAtEnd(removed1, "+M");

    //- Third: remove '>' if we have forward reaction
    //  and have product site
    if (site == "p" && !data.BR(data.nReac()))
    {
        removeFirstChar(tmp); 
    }

    word stochiometricFactor;

    word species;

    int startPos{0};

    int endPos{0};

    bool foundDigit{false};

    bool firstChar{false};

    bool extractSpecies{false};

    //- Reactant analyse
    for(unsigned int i=0; i<tmp.size(); i++)
    {
        if (!firstChar)
        {
            if
            (
                isdigit(tmp[i])
             || tmp[i] == '.'
            )
            {
                if (foundDigit)
                {
                    stochiometricFactor += tmp[i];
                }
                else
                {
                    stochiometricFactor = tmp[i];
                }

                foundDigit = true;
            }
            //- If first char is not a number, species name starts
            else
            {
                firstChar = true;

                startPos = i;
                
                if (!foundDigit)
                {
                    stochiometricFactor = "1";
                }
                else
                {
                    foundDigit = false;
                }
            }
        }

        {
            //- Move on till + comes or end of string reached
            if
            (
                tmp[i] == '+'
             || i == tmp.size()-1    
            )
            {
                if (i == tmp.size()-1)
                {
                    endPos = i+1;
                }
                else
                {
                    endPos = i;
                }

                firstChar = false;

                extractSpecies = true;
            }
        }

        //- Extract species
        if (extractSpecies)
        {
            species = tmp.substr(startPos, (endPos-startPos));

            extractSpecies = false;

            int nu{0};
            scalar tmp{0};

            if (site == "p")
            {
                tmp = stod(stochiometricFactor);

            }
            else if (site == "e")
            {
                tmp = stod(stochiometricFactor) * -1;
            }
            else
            {
                ErrorMsg
                (
                    "You only can call this function with 'e' or "
                    "'p' arguments.",
                    __FILE__,
                    __LINE__
                );
            }

            //- Check if value is integer, therefore we multiply the value
            //  by 10 (0.1 -> 1 | 1 -> 10) and make % operation
            const int modTmp = static_cast<int>(tmp * scalar(100));

            //- The modulo has to be zero if the input is an integer
            //  Input integer  -> 1   -> 10 % 10 = 0 
            //  Input integer  -> 2   -> 20 % 10 = 0
            //  Input !integer -> 1.2 -> 12 % 10 = 2
            //  Input !integer -> 0.4 -> 4  % 10 = 4
            if ((modTmp % 10) != 0)
            {
                ErrorMsg
                (
                    "The stochiometric factor is not an integer.\n    "
                    "Error occured in analyzing reaction site "
                    + reactionSite + " (" + site + ")",
                    __FILE__,
                    __LINE__
                );
            }
            else
            {
                //- Cast the scalar into an integer
                nu = static_cast<int>(tmp);
            }

            startPos = 0;
            endPos = 0;

            //- Store data
            data.speciesInReaction(species);

            if (site == "e")
            {
                data.nuEducts(species, nu);
                data.educt(species);
                data.forwardReactionOrder();
            }
            else if (site == "p")
            {
                data.nuProducts(species, nu);
                data.product(species);
                data.backwardReactionOrder();
            }
        }
    }
}


// ************************************************************************* //
