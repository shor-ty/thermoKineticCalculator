/*---------------------------------------------------------------------------*\
| =========                |                                                  |
| \\      /  A utomatic    | Holzmann-cfd                                     |
|  \\    /   F lamelet     | Version: 1.0                                     |
|   \\  /    C onstructor  | Web: www.Holzmann-cfd.de                         |
|    \\/                   |                                                  |
\*---------------------------------------------------------------------------*/
/*
»
»
»
»
»
»
»
»
»
»
\*---------------------------------------------------------------------------*/
//- system headers
#include <fstream>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <string>

//- user def. headers
#include "chemistry.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Chemistry::Chemistry()
:
    n_(-1),
    nDuplicate_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Chemistry::~Chemistry()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Chemistry::readChemkin
(
    const normalString& fileName
)
{
    //- read the whole file and store it
    stringField fileContent = openFile(fileName);

    //- loop through the fileContent
    forAll(fileContent, line)
    {
        //- STEP 1: remove all whitespaces
        stringField tmp = splitString(fileContent[line]);


        //- if not empty, and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            //- ELEMENTS BLOCK
            if
            (
                tmp[0] == "ELEMENTS"
             || tmp[0] == "ELEM"
            )
            {
                //- skip the line with the keyword ELEMENTS
                line++;

                //- loop till we reach the keyword END
                for (;;line++)
                {
                    tmp = splitString(fileContent[line]);

                    //- if line is not empty and no comment, proceed
                    if
                    (
                        !tmp.empty()
                     && tmp[0][0] != '!'
                    )
                    {
                        //- when the line content is END, leave that loop
                        if (tmp[0] == "END")
                        {
                            //- skip END
                            line++;
                            break;
                        }

                        //- when more elements are in one line
                        if (tmp.size() > 1)
                        {
                            forAll(tmp, element)
                            {
                                elements_.push_back(tmp[element]);
                            }
                        }
                        else
                        {
                            elements_.push_back(tmp[0]);
                        }
                    }
                }
            //- ELMENTS BLOCK END
            }


            //- SPECIES BLOCK
            else if
            (
                tmp[0] == "SPECIES"
            )
            {
                //- skip the line with the keyword SPECIES
                line++;

                //- loop till we reach the keyword END
                for (;;line++)
                {
                    tmp = splitString(fileContent[line]);

                    //- if line is not empty and no comment, proceed
                    if
                    (
                        !tmp.empty()
                     && tmp[0][0] != '!'
                    )
                    {
                        //- when the line content is END, leave that loop
                        if (tmp[0] == "END")
                        {
                            //- skip END
                            line++;
                            break;
                        }

                        //- when more SPECIES are in one line
                        if (tmp.size() > 1)
                        {
                            forAll(tmp, species)
                            {
                                species_.push_back(tmp[species]);
                            }
                        }
                        else
                        {
                            species_.push_back(tmp[0]);
                        }
                    }
                }
            //- SPECIES BLOCK END
            }


            //- THERMO BLOCK
            else if
            (
               (tmp[0] == "THERMO"
             && tmp[1] == "ALL")
             || tmp[0] == "THERMO"
            )
            {
                //- skip the line with the keyword THERMO
                line++;

                //- loop till we reach the keyword END
                for (;;line++)
                {
                    tmp = splitString(fileContent[line]);

                    //- if line is not empty and no comment, proceed
                    if
                    (
                        !tmp.empty()
                     && tmp[0][0] != '!'
                    )
                    {
                        //- when the line content is END, leave that loop
                        if (tmp[0] == "END")
                        {
                            //- skip END
                            line++;
                            break;
                        }
                    }
                }
            //- THERMO BLOCK END
            }


            //- REACTION BLOCK
            else if
            (
                tmp[0] == "REACTIONS"
            )
            {
                //- skip the line with the keyword THERMO
                line++;

                //- loop till we reach the keyword END
                for (;;line++)
                {
                    tmp = splitString(fileContent[line]);

                    //- if line is not empty and no comment, proceed
                    if
                    (
                        !tmp.empty()
                     && tmp[0][0] != '!'
                    )
                    {
                        //- when the line content is END, leave that loop
                        if (tmp[0] == "END")
                        {
                            //- skip END
                            line++;
                            break;
                        }

                        //- check for duplicate entrys
                        if (tmp[0] == "DUPLICATE")
                        {
                            //- skip next two lines (third with for loop)
                            line+= 2;

                            //- increment duplicate entry
                            nDuplicate_++;
                        }
                        else
                        {
                            //- check if '=' is in string (means reaction)
                            std::size_t found = fileContent[line].find('=');

                            if (found != std::string::npos)
                            {
                                //- increase reaction amount
                                n_++;

                                //- increment all matrixes
                                incrementMatrixesVectors();

                                //- save the elementar reaction
                                elementarReaction(fileContent[line]);

                                //- split elementar reaction into reactants
                                //  and product species get stochiometric
                                //  factors and update the matrix
                                //reactantsAndProducts();
                            }
                            else
                            std::cout << fileContent[line] << "\n";
                        }
                    }
                }
            //- REACTION BLOCK END
            }

        //- comment check
        }
    }
}



void Chemistry::update
(
    const normalString& str
)
{
    //- STEP 1: check if its a reaction
    std::size_t found = str.find("=");
    std::cout << str << "\n";

    //- STEP 2: if reaction found save:
    //
    //  +  a) elementar reaction
    //  +  b) arrhenius coeffs
    //  +  c) increment all matrixes and vectors
    //  +  d) increase the variable n_
    //  +  e) get stochiometric factors
    if (found != std::string::npos)
    {
        //- increment n_
        //  reason that we start with -1 is the first reaction
        //  --> becomes 0
        n_++;

    }


    //- STEP 1: save reaction
//
//
//    //- STEP 2: save arrhenius coeffs
//    arrheniusCoeffs(reaction);
//
//    //- STEP 3: check if THIRD BODY REACTION
//    thirdBodyReaction
//    (
//        fileContent,
//        line
//    );
//
//    //- STEP 4: handle THIRD BODY REACTION
//    if (TBR_[n_])
//    {
//        handleThirdBodyReaction
//        (
//            fileContent,
//            line
//        );
//    }
//
//    //- STEP 5: check if backward reaction is used
//    backwardReaction();
//

//
//    //- STEP 5: increment reaction counter
//    n_++;
}


void Chemistry::incrementMatrixesVectors()
{

    //- increment matrixes and vectors

        //- vector for saving reactions
        elementarReaction_.push_back("");

        //- matrix for stochiometric coeffs
        nu_.push_back(std::vector<double>(species_.size()));

        //- matrix of THIRD BODY M (composition of species)
        Mcomp_.push_back(std::vector<std::string>(0));

        //- matrix of THIRD BODY M (values of species)
        Mvalue_.push_back(std::vector<double>(0));

        //- matrix of ARRHENIUS coeffs
        arrheniusCoeffs_.push_back(std::vector<double>(3));

        //- matrix of TROE coeffs
        TROECoeffs_.push_back(std::vector<double>(5));

        //- matrix of ARRHENIUS coeffs for LOW pressure
        LOWCoeffs_.push_back(std::vector<double>(3));

        //- matrix of SRI coeffs
        SRICoeffs_.push_back(std::vector<double>(5));

        //- matrix of reactants that are used in reaction n_
        reactants_.push_back((std::vector<std::string>(3)));

        //- matrix of products that are used in reaction n_
        products_.push_back((std::vector<std::string>(3)));

        //- vector for THIRD BODY REACTION
        TBR_.push_back(false);

        //- vector for THIRD BODY REACTION LOW
        LOW_.push_back(false);

        //- vector for THIRD BODY REACTION TROE
        TROE_.push_back(false);

        //- vector for THIRD BODY REACTION SRI
        SRI_.push_back(false);

        //- vector for THIRD BODY REACTION of ENHANCEMENT FACTORS
        ENHANCE_.push_back(false);

        //- vector for backward reaction
        kb_.push_back(false);
}


void Chemistry::elementarReaction
(
    const normalString& reaction
)
{
    stringField tmp = splitString(reaction);
    normalString tmp2;

    //- re-arrange the string and remove arrhenius coeffs
    for (unsigned int i=0; i<tmp.size()-3; i++)
    {
        tmp2 += tmp[i];
    }

    elementarReaction_[n_] = tmp2;
}


void Chemistry::arrheniusCoeffs
(
    const normalString& reaction
)
{
    stringField tmp = splitString(reaction);

    //- not a complete check
    if (tmp.size() < 4)
    {
        std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
            << __LINE__ << " ++ Amount of entrys in "
            << elementarReaction_[n_] << " is bad. Found: " << tmp.size()
            << std::endl;
        std::terminate();
    }

    //- save arrhenius coeffs (last 3 entrys in tmp array)
    arrheniusCoeffs_[n_][0] = stod(tmp[tmp.size()-3]);
    arrheniusCoeffs_[n_][1] = stod(tmp[tmp.size()-2]);
    arrheniusCoeffs_[n_][2] = stod(tmp[tmp.size()-1]);
}


void Chemistry::thirdBodyReaction
(
    const stringField& fileContent,
    const unsigned int& line
)
{
    //- search for the indicator (only first match)

        //- for standard THIRD BODY REACTIONS
        //  definition:
        //  +  a) next line is a new reaction -> constant rate ind. of p
        //  +  b) next line is a enhancement factor
        normalString delimiter1 = "+M=";

        //- for modified THIRD BODY REACTIONS (TROE/LOW)
        normalString delimiter2 = "(+M)=";

        std::size_t found_standard = elementarReaction_[n_].find(delimiter1);
        std::size_t found_modified = elementarReaction_[n_].find(delimiter2);

    //- match means THIRD BODY REACTION
    //  definition:
    //  +  a) found_standard = true  and found_modified = false
    //  +  b) found_standard = false and found_modified = true
    //  +  c) both true  (not possible)
    //  +  d) both false (standard reaction)

    //- chase a)
    if
    (
        found_standard != std::string::npos
     && found_modified == std::string::npos
    )
    {
        TBR_[n_] = true;

        //- check next line for ENHANCE
        normalString delimiterE="/";
        std::size_t found_E = fileContent[line+1].find(delimiterE);

        //- match means, that ENHANCE factors found
        if (found_E != std::string::npos)
        {
            ENHANCE_[n_] = true;
        }

    }
    //- case b)
    else if
    (
        found_standard == std::string::npos
     && found_modified != std::string::npos
    )
    {
        TBR_[n_] = true;

        //- check next line for LOW
        normalString delimiterLOW="LOW";
        std::size_t found_LOW = fileContent[line+1].find(delimiterLOW);

        // match means, that LOW is found
        if (found_LOW != std::string::npos)
        {
            //- check next line for TROE
            normalString delimiterTROE="TROE";
            std::size_t found_TROE = fileContent[line+2].find(delimiterTROE);

            //- check next line for SRI
            normalString delimiterSRI="SRI";
            std::size_t found_SRI = fileContent[line+2].find(delimiterSRI);

            //- match means, that TROE or SRI was found
            if
            (
                found_TROE != std::string::npos
             || found_SRI != std::string::npos
            )
            {
                LOW_[n_] = true;

                //- TROE
                if (found_TROE != std::string::npos)
                {
                    TROE_[n_] = true;
                }
                else
                {
                    SRI_[n_] = true;
                }

                //- check next line for ENHANCE
                normalString delimiterE="/";
                std::size_t found_E = fileContent[line+2].find(delimiterE);

                //- match means, that ENHANCE factors found
                if (found_E != std::string::npos)
                {
                    ENHANCE_[n_] = true;
                }
            }
            //- only LOW is used
            else
            {
                LOW_[n_] = true;

                //- check next line for ENHANCE
                normalString delimiterE="/";
                std::size_t found_E = fileContent[line+2].find(delimiterE);

                //- match means, that ENHANCE factors found
                if (found_E != std::string::npos)
                {
                    ENHANCE_[n_] = true;
                }
            }
        }
        //- something went wrong
        else
        {
            //- not possible ERROR
            std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
                << __LINE__ << " ++ THIRD BODY REACTION | LOW not found\n"
                << " ++ REACTION = " << elementarReaction_[n_] << std::endl;
            std::terminate();
        }
    }
    //- case c)
    else if
    (
        found_standard != std::string::npos
     && found_modified != std::string::npos
    )
    {
        //- not possible ERROR
        std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
            << __LINE__ << " ++ THIRD BODY REACTION problem" << std::endl;
        std::terminate();
    }
}


void Chemistry::handleThirdBodyReaction
(
    const stringField& fileContent,
    const unsigned int& line
)
{
    //- EXCEPTION 1
    //  only ENHANCED factors
    if
    (
        ENHANCE_[n_]
     && !LOW_[n_]
     && !TROE_[n_]
     && !SRI_[n_]
    )
    {
        //- update the ENHANCED factors matrix
        enhancedFactors(fileContent[line+1]);

    }

    //- EXCEPTION 2
    //  only LOW and ENHANCED factors
    else if
    (
        ENHANCE_[n_]
     && LOW_[n_]
     && !TROE_[n_]
     && !SRI_[n_]
    )
    {
        //- update the LOW pressure arrhenius coeffs matrix
        LOWCoeffs(fileContent[line+1]);

        //- update the ENHANCED factors matrix
        enhancedFactors(fileContent[line+2]);
    }

    //- EXCEPTION 3
    //  only LOW TROE AND ENHANCED
    else if
    (
        ENHANCE_[n_]
     && LOW_[n_]
     && TROE_[n_]
     && !SRI_[n_]
    )
    {
        //- update the LOW pressure arrhenius coeffs matrix
        LOWCoeffs(fileContent[line+1]);

        //- update the TROE coeffs matrix
        TROECoeffs(fileContent[line+2]);

        //- update the ENHANCED factors matrix
        enhancedFactors(fileContent[line+3]);
    }
    //- EXCEPTION 4
    //  only LOW SRI AND ENHANCED
    else if
    (
        ENHANCE_[n_]
     && LOW_[n_]
     && !TROE_[n_]
     && SRI_[n_]
    )
    {
        //- update the LOW pressure arrhenius coeffs matrix
        LOWCoeffs(fileContent[line+1]);

        //- update the TROE coeffs matrix
        SRICoeffs(fileContent[line+2]);

        //- update the ENHANCED factors matrix
        enhancedFactors(fileContent[line+3]);
    }
    //- something went wrong
    //  TODO here no LOW, TROE, SRI, ENHANCED factors for THREE BODY
    else
    {
        std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
            << __LINE__ << " ++ Something is wrong with reaction: "
            << elementarReaction_[n_] << std::endl;
        std::terminate();
    }
}


void Chemistry::enhancedFactors
(
    const normalString& enhancedFactors
)
{
    //- remove all white space and split
    stringField tmp = splitString(enhancedFactors);

    normalString enhanced;

    forAll(tmp, i)
    {
        enhanced += tmp[i];
    }

    //- split enhanced string with delimiter '/'
    tmp = splitString(enhanced, '/');

    forAll(tmp, i)
    {
        //- modul
        unsigned int modul=i%2;

        //- value
        if (modul)
        {
            Mvalue_[n_].push_back(stod(tmp[i]));
        }
        //- species
        else
        {
            Mcomp_[n_].push_back(tmp[i]);
        }
    }
}


void Chemistry::LOWCoeffs
(
    const normalString& LOWcoeffs
)
{
    //- find first '/'
    normalString delimiter="/";
    std::size_t found = LOWcoeffs.find(delimiter);

    //- remove all letters from 0 till pos of '/'
    normalString tmp = LOWcoeffs.substr(found+1,LOWcoeffs.size());

    //- find second '/'
    found = tmp.find(delimiter);

    //- remove all letters from pos of '/' till end
    normalString tmp2 = tmp.substr(0,found);

    //- remove all white space and split
    stringField coeffs = splitString(tmp2);

    //- check if 3 values are available
    if (coeffs.size() != 3)
    {
        std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
            << __LINE__ << " ++ More or less than 3 arrhenius coeffs in "
            << elementarReaction_[n_] << ". Found: " << coeffs.size()
            << std::endl;
        std::terminate();
    }

    //- update matrix
    forAll(coeffs, i)
    {
        LOWCoeffs_[n_][i] = stod(coeffs[i]);
    }
}


void Chemistry::TROECoeffs
(
    const normalString& TROEcoeffs
)
{
    //- find first '/'
    normalString delimiter="/";
    std::size_t found = TROEcoeffs.find(delimiter);

    //- remove all letters from 0 till pos of '/'
    normalString tmp = TROEcoeffs.substr(found+1,TROEcoeffs.size());

    //- find second '/'
    found = tmp.find(delimiter);

    //- remove all letters from pos of '/' till end
    normalString tmp2 = tmp.substr(0,found);

    //- remove all white space and split
    stringField coeffs = splitString(tmp2);

    //- check if more than values are available
    if (coeffs.size() > 5)
    {
        std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
            << __LINE__ << " ++ More than 4 TROE coeffs in "
            << elementarReaction_[n_] << ". Found: " << coeffs.size()
            << std::endl;
        std::terminate();
    }

    //- update matrix
    forAll(coeffs, i)
    {
        if (!coeffs[i].empty())
        {
            TROECoeffs_[n_][i] = stod(coeffs[i]);
        }
    }
}


void Chemistry::SRICoeffs
(
    const normalString& SRIcoeffs
)
{
    //- find first '/'
    normalString delimiter="/";
    std::size_t found = SRIcoeffs.find(delimiter);

    //- remove all letters from 0 till pos of '/'
    normalString tmp = SRIcoeffs.substr(found+1,SRIcoeffs.size());

    //- find second '/'
    found = tmp.find(delimiter);

    //- remove all letters from pos of '/' till end
    normalString tmp2 = tmp.substr(0,found);

    //- remove all white space and split
    stringField coeffs = splitString(tmp2);

    //- check if more than values are available
    if (coeffs.size() > 5)
    {
        std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
            << __LINE__ << " ++ More than 5 SRI coeffs in "
            << elementarReaction_[n_] << ". Found: " << coeffs.size()
            << std::endl;
        std::terminate();
    }

    //- update matrix
    forAll(coeffs, i)
    {
        if (!coeffs[i].empty())
        {
            SRICoeffs_[n_][i] = stod(coeffs[i]);
        }
    }
}


void Chemistry::backwardReaction()
{
    normalString delimiter1 = "=";
    normalString delimiter2 = "<";
    normalString delimiter3 = ">";

    //- reactions definition
    //
    //  +  a)  A+A=B   (forward and backward)
    //  +  b)  A+A<=>B (forward and backward)
    //  +  c)  A+A=>B  (only forward reaction)
    //  +  d)  A+A<=B  (only backward reaction)
    //  d) is not implemented

    std::size_t found1 = elementarReaction_[n_].find(delimiter1);
    std::size_t found2 = elementarReaction_[n_].find(delimiter2);
    std::size_t found3 = elementarReaction_[n_].find(delimiter3);

    //- handle only stuff which include backward reaction because
    //  standard setting is false

    //- a)
    if
    (
        found1 != std::string::npos
     && found2 == std::string::npos
     && found3 == std::string::npos
    )
    {
        kb_[n_] = true;
    }
    //- b)
    else if
    (
        found1 != std::string::npos
     && found2 != std::string::npos
     && found3 != std::string::npos
    )
    {
        kb_[n_] = true;
    }
    //- if for some reason other signs than
    //  =
    //  =>
    //  <=>
    //  are used, programm will canceled in fucntion
    //  reactantsAndProducts()
}


void Chemistry::reactantsAndProducts()
{
    //- STEP 1: split into reactants and products
    normalString delimiter1 = "=";
    normalString delimiter2 = ">";
    normalString delimiter3 = "<";
    normalString delimiter4 = "";

    //- reactions definition
    //
    //  +  a)  A+A=B   (forward and backward)
    //  +  b)  A+A<=>B (forward and backward)
    //  +  c)  A+A=>B  (only forward reaction)
    //  +  d)  A+A<=B  (only backward reaction)
    //  d) is not implemented
    {
        std::size_t found1 = elementarReaction_[n_].find(delimiter1);
        std::size_t found2 = elementarReaction_[n_].find(delimiter2);
        std::size_t found3 = elementarReaction_[n_].find(delimiter3);

        //- handle only stuff which include backward reaction because
        //  standard setting is false

        //- a)
        if
        (
            found1 != std::string::npos
         && found2 == std::string::npos
         && found3 == std::string::npos
        )
        {
            delimiter4 = "=";
        }
        //- b)
        else if
        (
            found1 != std::string::npos
         && found2 != std::string::npos
         && found3 != std::string::npos
        )
        {
            delimiter4 = "<=>";
        }
        //- c)
        else if
        (
            found1 != std::string::npos
         && found2 != std::string::npos
         && found3 == std::string::npos
        )
        {
            delimiter4 = "=>";
        }
        else
        {
            std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
                << __LINE__ << " ++ Reaction direction is not clear for "
                << elementarReaction_[n_] << std::endl;
            std::terminate();
        }
    }

    normalString reactants =
        elementarReaction_[n_].substr
        (
            0,
            elementarReaction_[n_].find(delimiter4)
        );

    normalString products =
        elementarReaction_[n_].substr
        (
            elementarReaction_[n_].find(delimiter4)+delimiter4.size(),
            elementarReaction_[n_].size()
        );

    //- STEP 2: remove THIRD BODY
    if (TBR_[n_])
    {
        removeThirdBody(reactants);
        removeThirdBody(products);
    }

    //- STEP 3: split reactants and products into species
    stringField reacElements = splitString(reactants, '+');
    stringField prodElements = splitString(products, '+');

    //- STEP 4: update stochiometric
    updateStochiometricMatrix(reacElements, 0);
    updateStochiometricMatrix(prodElements, 1);
}


void Chemistry::updateStochiometricMatrix
(
    const stringField& speciesField,
    const unsigned int x
)
{
    //- bool variable for match information
    bool found{false};


    //- loop through the field
    forAll(speciesField, species)
    {
        found = false;

        //- loop through all single letter of the species
        for (unsigned int pos = 0; pos < speciesField[species].size(); pos++)
        {
            //- compare single letter with ASCII table (LETTERS)
            for (unsigned int j = 65; j <= 90; j++)
            {
                char c = static_cast<char>(j);

                if (c == speciesField[species][pos])
                {
                    found = true;
                    break;
                }
            }

            //- if letter found, check which position
            if (found)
            {
                //- pos = 0 (first letter) --> nu = 1
                if (pos == 0)
                {
                    updateStochiometricMatrix
                    (
                        speciesField[species],
                        1,
                        x
                    );
                }
                //- if not pos = 0, scalar has to be splitted
                else
                {
                    //- line below also work
                    //  scalar nu = stod(speciesField[species]);
                    scalar nu = stod(speciesField[species].substr(0,pos));
                    normalString specie =
                        speciesField[species].substr
                        (
                            pos,speciesField[species].size()
                        );

                    updateStochiometricMatrix
                    (
                        specie,
                        nu,
                        x
                    );
                }

                //- move to next species
                break;
            }
        //- loop over all single letter in species
        }
    //- loop over all species
    }
}


void Chemistry::updateStochiometricMatrix
(
    const normalString& species,
    const scalar nu,
    const unsigned int x
)
{
    //- STEP 1: find id of species
    unsigned int ID{0};
    bool found{false};

    forAll(species_, id)
    {
        if (species == species_[id])
        {
            found = true;
            ID = id;
        }
    }

    //- check if species found or not
    if (!found)
    {
            std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
                << __LINE__ << " ++ Species " << species << " not defined in "
                << "the SPECIES section of the chemkin file." << std::endl;
            std::terminate();
    }


    //- STEP 2: get actual nu out from the matrix
    scalar nuTmp = nu_[n_][ID];

    //- STEP 3: check if product or reactant and add the new value
    //  needed because of HO+HO (twice)

        //- reactants
        if (x == 0)
        {
            nuTmp -= nu;
        }
        //- products
        else if (x == 1)
        {
            nuTmp += nu;
        }

    //- STEP 4: save the new stochiometric factor
    nu_[n_][ID] = nuTmp;
}


void Chemistry::removeThirdBody
(
    normalString& reaction
)
{
    //- delimiter
    normalString delimiter = "";

    //- check which THIRD BODY reaction is used
    //  TBR true
    //  && (LOW true || TROE ture || SRI true)
    if
    (
        TBR_[n_]
     && (LOW_[n_] || TROE_[n_] || SRI_[n_])
    )
    {
        //- THIRD BODY KEYWORD (+M)
        delimiter = "(+M)";

    }
    //  TBR true
    //  && LOW false && TROE false && SRI false
    else if
    (
        TBR_[n_]
     && !LOW_[n_]
     && !TROE_[n_]
     && !SRI_[n_]
    )
    {
        //- THIRD BODY KEYWORD +M
        delimiter = "+M";
    }

    //- find the match
    std::size_t found_TB = reaction.find(delimiter);

    //- modify reaction (remove TB)
    reaction = reaction.substr(0,found_TB);
}


stringField Chemistry::splitString
(
    const normalString& str
)
{
    //- split string, delimiter is whitespace
    std::istringstream tmp(str);

    stringField strArray_
    {
        std::istream_iterator<std::string>{tmp},
        std::istream_iterator<std::string>{}
    };

    return strArray_;
}


stringField Chemistry::splitString
(
    const normalString& str,
    const char delimiter
)
{
    //- split the line with delimiter
    std::stringstream tmp(str);
    normalString element;
    stringField elements;

    while (std::getline(tmp, element, delimiter))
    {
        elements.push_back(element);
    }

    //- return the field
    return elements;
}


stringField Chemistry::openFile
(
    const normalString& fileName
)
{
    //- new object
    std::ifstream file;

    //- open file
    file.open(fileName.c_str(), std::ios::in);

    //- check file
    if (!file.good())
    {
        std::cerr<< "\n ++ ERROR: Could not open file \""
                 << fileName << "\"..."
                 << "\n ++ Error occur in file " << __FILE__
                 << " line " << __LINE__ << std::endl;
        std::terminate();
    }

    //- stored line
    normalString fileLine;

    //- file content
    stringField fileContent(0);

    //- read the whole file and store the lines inside fileContent
    while (!file.eof())
    {
        std::getline(file, fileLine);
        fileContent.push_back(fileLine);
    }

    return fileContent;
}


void Chemistry::summary() const
{

    //- THIRD BODY REACTIONS
    unsigned int TBR{0};
    unsigned int LOW{0};
    unsigned int TROE{0};
    unsigned int SRI{0};
    unsigned int ENH{0};

    forAll(TBR_, i)
    {
        if (TBR_[i])
        {
            TBR++;
        }
        if (LOW_[i])
        {
            LOW++;
        }
        if (TROE_[i])
        {
            TROE++;
        }
        if (SRI_[i])
        {
            SRI++;
        }
        if (ENHANCE_[i])
        {
            ENH++;
        }
    }

    std::cout<< "----------------------------------------------------------\n"
             << "                   CHEMISTRY SUMMARY                      \n"
             << "----------------------------------------------------------\n"
             << std::left
             << std::setw(40) << " No. of elements: "
             << std::setw(10) << elements_.size() << "\n"
             << std::setw(40) << " No. of species in reactions: "
             << std::setw(20) << species_.size() << "\n"
             << std::setw(40) << " No. of elementar reactions: "
             << std::setw(20) << n_+1 << "\n"
             << std::setw(40) << " No. of third body reactions: "
             << std::setw(20) << TBR << "\n"
             << std::setw(40) << " No. of low pressure reactions: "
             << std::setw(20) << LOW-TROE-SRI << "\n"
             << std::setw(40) << " No. of TROE reactions: "
             << std::setw(20) << TROE << "\n"
             << std::setw(40) << " No. of SRI reactions: "
             << std::setw(20) << SRI << "\n"
             << std::setw(40) << " No. of enhanced factor reactions: "
             << std::setw(20) << ENH-LOW << "\n"
             << std::setw(40) << " No. of dublicated reactions: "
             << std::setw(20) << nDuplicate_ << "\n"
             << "----------------------------------------------------------\n";
}


void Chemistry::createReactionRateMatrix()
{
    //- loop through all species
    forAll(species_, id)
    {
        reactionRateMatrix_.push_back(std::vector<double>(0));

        //- loop through the nu matrix to figure out which elementar reaction
        //  influences the species
        forAll(nu_, r)
        {
            if (nu_[r][id] != 0)
            {
                reactionRateMatrix_[id].push_back(r);
            }
        }
    }

}
