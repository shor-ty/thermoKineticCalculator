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
    nDuplicate_(0),
    themodynamic_(false)
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
                            break;
                        }

                        //- when more SPECIES are in one line
                        if (tmp.size() > 1)
                        {

                            forAll(tmp, species)
                            {
                                //- put species at the end of the vector
                                species_.push_back(tmp[species]);

                                //- increment all thermodynamic variables
                                //  matrixes, vectors etc.
                                thermodynamicDataIncrement();

                                //- increment all transport vectors
                                transportDataIncrement();
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
                //- set the variable
                themodynamic_ = true;
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

                                //- backward reaction
                                backwardReaction();

                                //- check if THIRD BODY REACTION
                                found = fileContent[line].find("+M");

                                if (found != std::string::npos)
                                {
                                    TBR_[n_] = true;
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
                                    LOWCoeffs(fileContent[line]);
                                    LOW_[n_] = true;
                                }
                                //- TROE parameters
                                else if (foundTROE != std::string::npos)
                                {
                                    TROECoeffs(fileContent[line]);
                                    TROE_[n_] = true;
                                }
                                //- SRI parameters
                                else if (foundSRI != std::string::npos)
                                {
                                    SRICoeffs(fileContent[line]);
                                    SRI_[n_] = true;
                                }
                                else
                                {
                                    enhancedFactors(fileContent[line]);
                                    ENHANCE_[n_] = true;
                                }
                            }
                        //- working else
                        }
                    //- check for empty and comments
                    }
                //- loop in reaction block
                }
            //- REACTION BLOCK END
            }
        //- comment check
        }
    //- loop through the file
    }


    //- split elementar reaction into reactants and product
    //  get stochiometric factors and update the matrix
    reactantsAndProducts();
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
    forAll(elementarReaction_, i)
    {
        std::size_t found1 = elementarReaction_[i].find(delimiter1);
        std::size_t found2 = elementarReaction_[i].find(delimiter2);
        std::size_t found3 = elementarReaction_[i].find(delimiter3);

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
                << elementarReaction_[i] << std::endl;
            std::terminate();
        }

        normalString reactants =
            elementarReaction_[i].substr
            (
                0,
                elementarReaction_[i].find(delimiter4)
            );

        normalString products =
            elementarReaction_[i].substr
            (
                elementarReaction_[i].find(delimiter4)+delimiter4.size(),
                elementarReaction_[i].size()
            );

        //- STEP 2: remove THIRD BODY
        if (TBR_[i])
        {
            removeThirdBody(reactants, i);
            removeThirdBody(products, i);
        }

        //- STEP 3: split reactants and products into species
        stringField reacElements = splitString(reactants, '+');
        stringField prodElements = splitString(products, '+');

        //- STEP 4: update stochiometric
        updateStochiometricMatrix(reacElements, 0, i);
        updateStochiometricMatrix(prodElements, 1, i);
    }
}


void Chemistry::updateStochiometricMatrix
(
    const stringField& speciesField,
    const unsigned int x,
    const unsigned int& i
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
                        x,
                        i
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
                        x,
                        i
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
    const unsigned int x,
    const unsigned int& i
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
    scalar nuTmp = nu_[i][ID];

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
    nu_[i][ID] = nuTmp;
}


void Chemistry::removeThirdBody
(
    normalString& reaction,
    const unsigned int& i
)
{
    //- delimiter
    normalString delimiter = "";

    //- check which THIRD BODY reaction is used
    //  TBR true
    //  && (LOW true || TROE ture || SRI true)
    if
    (
        TBR_[i]
     && (LOW_[i] || TROE_[i] || SRI_[i])
    )
    {
        //- THIRD BODY KEYWORD (+M)
        delimiter = "(+M)";

    }
    //  TBR true
    //  && LOW false && TROE false && SRI false
    else if
    (
        TBR_[i]
     && !LOW_[i]
     && !TROE_[i]
     && !SRI_[i]
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
             << std::setw(20) << LOW << "\n"
             << std::setw(40) << " No. of TROE reactions: "
             << std::setw(20) << TROE << "\n"
             << std::setw(40) << " No. of SRI reactions: "
             << std::setw(20) << SRI << "\n"
             << std::setw(40) << " No. of enhanced factor reactions: "
             << std::setw(20) << ENH << "\n"
             << std::setw(40) << " No. of dublicated reactions: "
             << std::setw(20) << nDuplicate_ << "\n"
             << "----------------------------------------------------------\n";


//             forAll(elementarReaction_, i)
//             {
//                std::cout << elementarReaction_[i] << "  >> ";
//                forAll(nu_[i], a)
//                {
//                    std::cout <<  nu_[i][a] << "   " ;
//                }
//                std::cout << "\n";
//             }
}

void Chemistry::readChemKinThermo
(
    const normalString& fileName
)
{
    //- output
    std::cout << "Reading thermodynamic properties (" << fileName << ")\n";

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

                        readNASA
                        (
                            fileContent,
                            line,
                            species_
                        );
                    }
                }
            }
        }
    }
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


void Chemistry::checkThermo
(
    const normalString& fileName
) const
{
    forAll(species_, i)
    {
        if (!NASA(i))
        {
            std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
                << __LINE__ << " ++ Species " << species_[i] << " has no "
                << "thermodynamic. No entry found in thermodynamic file "
                << fileName << std::endl;
            std::terminate();
        }
    }
    std::cout<< " ++ Thermodynamic database O.K.\n\n";
}

void Chemistry::checkTrans
(
    const normalString& fileName
) const
{
    forAll(species_, i)
    {
        if (!TRANS(i))
        {
            std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
                << __LINE__ << " ++ Species " << species_[i] << " has no "
                << "transport properties. No entry found in transport file "
                << fileName << std::endl;
            std::terminate();
        }
    }
    std::cout<< " ++ Transport database O.K.\n\n";
}


bool Chemistry::thermo() const
{
    return themodynamic_;
}


stringField Chemistry::species() const
{
    return species_;
}
