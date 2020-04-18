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
    Info<< " c-o Reading chemistry data\n"
        << "     >> " << file_ << "\n" << endl;

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
    if (debug_)
    {
        Info<< "Start reading the ELEMENT block\n" << endl;
    }

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

    if (debug_)
    {
        Info<< "End reading the ELEMENT block\n" << endl;
    }
}


void AFC::ChemistryReader::readSpeciesBlock
(
   const stringList& fileContent,
   ChemistryData& data
)
{
    if (debug_)
    {
        Info<< "Start reading the SPECIES block\n" << endl;
    }

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

    if (debug_)
    {
        Info<< "End reading the SPECIES block\n" << endl;
    }
}


void AFC::ChemistryReader::readThermoBlock
(
    const stringList& fileContent,
    ChemistryData& data
)
{
    if (debug_)
    {
        Info<< "Start reading the THERMO block\n" << endl;
    }

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
        // TODO
        //data.setThermo();
    }

    if (debug_)
    {
        Info<< "End reading the THERMO block\n" << endl;
    }
}

void AFC::ChemistryReader::readReactionBlock
(
    const stringList& fileContent,
    ChemistryData& data
)
{
    if (debug_)
    {
        Info<< "Start reading the REACTION block\n" << endl;
    }

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

    //- For duplicated entries
    bool duplicate{false};

    //- Reading REACTION block
    for (unsigned int line = lineNoKeyword+1; line < lineNoEnd; line++)
    {
        string tString = fileContent[line];

        //- Remove all comments (!)
        removeComment(tString);
        stringList tmp = splitStrAtWS(tString);

        if (!tmp.empty())
        {
            //- Check for duplicate entrys (not taken into consideration)
            if (tmp[0] == "DUPLICATE" || tmp[0] == "DUP")
            {
                //- Found first duplicated key
                if (!duplicate)
                {
                    duplicate = true;

                    //- Increment duplicate entry
                    data.incrementDuplicated();
                }
                else
                {
                    duplicate = false;
                }
            }
            else if (!duplicate)
            {
                //- Check if '=' is in string (means reaction)
                std::size_t found = fileContent[line].find('=');

                if (found != std::string::npos)
                {
                    data.incrementReac();
                    data.incrementMatrixesVectors();

                    //- If collision partner found >> in backets (+M), (+H2)
                    //  This donates a third body reaction (TBR)
                    std::size_t found1 = tmp[0].find('(');
                    std::size_t found2 = tmp[0].find(')');
                    std::size_t found3 = fileContent[line].find("+M");

                    if
                    (
                        (found1 != std::string::npos)
                     && (found2 != std::string::npos)
                    )
                    {
                        //- Set reaction to be a TBR
                        data.TBR(true);

                        //- Extract collision partner
                        data.collisionPartner
                        (
                            extractBetweenKeys(fileContent[line])
                        );
                    }
                    else if (found3 != std::string::npos)
                    {
                        //- Set reaction to be a TBR
                        data.TBR(true);
                        data.collisionPartner("+M");
                    }

                    //- Analyze the reaction
                    analyzeReaction(tmp, data);
                }
                // No reaction, something different
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

    if (debug_)
    {
        Info<< "End reading the REACTION block\n" << endl;
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
            if (!tmp.empty() && tmp[0] == i)
            {
                start = lineNo;
                break;
            }
        }

        //- Search end after keyword found
        if (!tmp.empty() && tmp[0] == "END" && start != -1)
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
    const stringList& line,
    ChemistryData& data
)
{
    const string& reaction = line[0];
    data.elementarReaction(reaction);

    //- STEP 1: analyze reaction
    string delimiter1 = "=";
    string delimiter2 = "<";
    string delimiter3 = ">";

    //- Reactions definition
    //
    //  +  a)  A+A=B   (forward and backward)
    //  +  b)  A+A<=>B (forward and backward)
    //  +  c)  A+A=>B  (only forward reaction)
    //  +  d)  A+A<=B  (only backward reaction)

    size_t found1 = reaction.find(delimiter1);
    size_t found2 = reaction.find(delimiter2);
    size_t found3 = reaction.find(delimiter3);

    //- a)
    if
    (
        found1 != std::string::npos
     && found2 == std::string::npos
     && found3 == std::string::npos
    )
    {
        data.FR(true);
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
        data.FR(true);
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
        data.FR(true);
        data.BR(false);
    }
    else if
    (
        found1 != std::string::npos
     && found2 != std::string::npos
     && found3 == std::string::npos
    )
    {
        data.FR(false);
        data.BR(true);
    }

    // STEP 2: insert arrhenius coeffs
    data.arrheniusCoeffs
    (
        stod(line[1]),
        stod(line[2]),
        stod(line[3])
    );

    // STEP 3: get stochiometric values
    //  + educts negativ
    //  + products positiv

    // a) split into educts and products
    const stringList sites = splitStrAtDelimiter(reaction, '=');

    const string educ = sites[0];
    const string prod = sites[1];

    // b) Analyze the educt and product site according to the species
    //    and its stochiometric numbers
    {
        //- Educt site
        analyzeReacSite(educ, "e", data);

        //- Product site
        analyzeReacSite(prod, "p", data);
    }

    //- All stochiometric coefficients for species stored. Hence we
    //  are able to calculate the global reaction order and exponents
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
    string reactionSite,
    const word site,
    ChemistryData& data
)
{
    //- First: check if actual reaction is a TBR, if yes remove the
    //  third body reaction partner

    if (data.TBR(data.nReac()))
    {
        const word collPartner = data.collisionPartner(data.nReac());

        //- Remove collision partner from reaction (partner is at the end)
        //  TODO check if collision partner is not at the end to be more dynamic
        //  However, not sure if the convention allows to set it somewhere else
        reactionSite = removeAtEnd(reactionSite, collPartner);
    }

    //- Second: check if '<' or '>' is inside the reaction
    //  a) for educt site we can have '<'
    //  b) for product site we can have '>'
    std::size_t found{false};

    if (site == "e")
    {
        found = reactionSite.find('<');

        if (found != std::string::npos)
        {
            removeLastChar(reactionSite);
        }
    }
    else if (site == "p")
    {
        found = reactionSite.find('>');

        if (found != std::string::npos)
        {
            removeFirstChar(reactionSite);
        }
    }
    else
    {
        ErrorMsg
        (
            "    You only can call this function with 'e' or "
            "'p' arguments.",
            __FILE__,
            __LINE__
        );
    }

    word stochiometricFactor;

    word species;

    int startPos{0};

    int endPos{0};

    bool foundDigit{false};

    bool firstChar{false};

    bool extractSpecies{false};

    //- Reactant analyse
    for(unsigned int i=0; i<reactionSite.size(); i++)
    {
        if (!firstChar)
        {
            //- Already prepared for non-int (more complex) reactions
            //  TODO complete implementation
            if (isdigit(reactionSite[i]) || reactionSite[i] == '.')
            {
                if (foundDigit)
                {
                    stochiometricFactor += reactionSite[i];
                }
                else
                {
                    stochiometricFactor = reactionSite[i];
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
            if (reactionSite[i] == '+' || i == reactionSite.size()-1)
            {
                if (i == reactionSite.size()-1)
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
            species = reactionSite.substr(startPos, (endPos-startPos));

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

            //- Check if value is an integer, therefore we multiply the value
            //  by 10 (0.1 -> 1 | 1 -> 10) and make % operation
            const int modTmp = static_cast<int>(tmp * scalar(10));

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
