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
    n_(0),
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
        //- ELEMENTS
        //  get the elements, that are included
        if (fileContent[line] == "ELEMENTS")
        {
            //- skip the line with the keyword ELEMENTS
            line++;

            //- loop till we reach the keyword END
            for (;;line++)
            {
                //- EXEPTION 1
                //  when the line content is empty
                if (fileContent[line].empty())
                {
                    //- skip that line
                    line++;
                }

                //- EXEPTION 2
                //  when the line content is a comment »!«
                stringField lineArray = splitString(fileContent[line]);

                //- when no comment, proceed
                if (lineArray[0] != "!")
                {
                    //- when the line content is END, leave that loop
                    if (fileContent[line] == "END")
                    {
                        //- skip END
                        line++;
                        break;
                    }

                    //- EXCEPTION 3
                    //  when more elements are in one line
                    if (lineArray.size() > 1)
                    {
                        forAll(lineArray, element)
                        {
                            elements_.push_back(lineArray[element]);
                        }
                    }
                    else
                    {
                        elements_.push_back(lineArray[0]);
                    }
                }
            }
        //- ELEMENTS end
        }


        //- SPECIES
        //  get the species, that are included
        if (fileContent[line] == "SPECIES")
        {
            //- skip the line with the keyword SPECIES
            line++;

            //- loop till we reach the keyword END
            for (;;line++)
            {
                //- EXEPTION 1
                //  when the line content is empty
                if (fileContent[line].empty())
                {
                    //- skip that line
                    line++;
                }

                //- EXCEPTION 2
                //  when the line content is a comment »!«
                stringField lineArray = splitString(fileContent[line]);

                //- when no comment, proceed
                if (lineArray[0] != "!")
                {
                    //- when the line content is END, leave that loop
                    if (fileContent[line] == "END")
                    {
                        //- skip END
                        line++;
                        break;
                    }

                    //- EXCEPTION 3
                    //  when more SPECIES are in one line
                    if (lineArray.size() > 1)
                    {
                        forAll(lineArray, species)
                        {
                            species_.push_back(lineArray[species]);
                        }
                    }
                    else
                    {
                        species_.push_back(lineArray[0]);
                    }
                }
            }
        //- SPECIES end
        }


        //- REACTIONS
        //  get the reactions that are included
        if (fileContent[line] == "REACTIONS")
        {
            //- skip the line with the keyword ELEMENTS
            line++;

            //- variables used in that section
            std::size_t found;

            //- loop till we reach the keyword END
            for (;;line++)
            {
                //- EXCEPTION 1
                //  when the line content is empty
                if (fileContent[line].empty())
                {
                    //- skip that line
                    line++;
                }

                //- EXCEPTION 2
                //  when the line content is a comment »!«
                stringField lineArray = splitString(fileContent[line]);

                //- EXCEPTION 3
                //  when the line content contains 'DUPLICATE'
                //  skip the next 3 lines
                //  TODO - get the point!
                if (lineArray[0] == "DUPLICATE")
                {
                    //- skip 3 lines
                    line+=3;

                    //- increment duplicated entrys
                    nDuplicate_++;

                    //- re-initialize the array with the new line
                    lineArray = splitString(fileContent[line]);
                }

                //- when no comment, proceed
                if (lineArray[0] != "!")
                {
                    //- when the line content is END, leave that loop
                    if (fileContent[line] == "END")
                    {
                        //- skip END
                        line++;
                        break;
                    }

                    //- check if the line contains '=' sign
                    found = fileContent[line].find('=');

                    //- the line content is a reaction
                    if (found != std::string::npos)
                    {
                        //- tmp
                        normalString reaction = fileContent[line];

                        //- increment the size of all vectors and matrixes
                        incrementMatrixesVectors();

                        //- main functionallity
                        //
                        //  +  1: take the reaction string
                        //  +  2: handle exceptions (TROE/LOW)
                        //  +  3: get species and stochiometric factors
                        update
                        (
                            reaction,
                            fileContent,
                            line
                        );
                    }
                }
            }
        //- REACTIONS end
        }
    //- loop through the file content end
    }
}



void Chemistry::update
(
    const normalString& reaction,
    const stringField& fileContent,
    const unsigned int& line
)
{

    //- STEP 1: save reaction
    elementarReaction(reaction);

    //- STEP 2: arrhenius coeffs
    arrheniusCoeffs(reaction);

    //- STEP 3: check if THIRD BODY REACTION
    thirdBodyReaction
    (
        fileContent,
        line
    );

    //- s

    //- increment reaction counter
    n_++;
}


void Chemistry::incrementMatrixesVectors()
{

    //- increment matrixes and vectors

        //- vector for saving reactions
        elementarReaction_.push_back("");

        //- matrix for stochiometric coeffs
        nu_.push_back(std::vector<double>(species_.size()));

        //- matrix of THIRD BODY M (composition of species)
        Mcomp_.push_back(std::vector<std::string>(species_.size()));

        //- matrix of THIRD BODY M (values of species)
        Mvalue_.push_back(std::vector<double>(species_.size()));

        //- matrix of ARRHENIUS coeffs
        arrheniusCoeffs_.push_back(std::vector<double>(3));

        //- matrix of TROE coeffs
        TROECoeffs_.push_back(std::vector<double>(3));

        //- matrix of ARRHENIUS coeffs for LOW pressure
        LOWCoeffs_.push_back(std::vector<double>(3));

        //- vector of fall off reactions (TROE) (0:=no | 1:=yes)
        nTROE_.push_back(0);

        //- vector of low pressure reactions
        nLOW_.push_back(0);

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

    //- check if tmp array is minimum size 4
    if (tmp.size() < 4)
    {
        std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
            << __LINE__ << " ++ Array size is less than 4" << std::endl;
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


normalString Chemistry::splitReaction
(
    const normalString& str
)
{
    int i = 0;
    //- split string, delimiter is:
    //  =
    //  =>
    //  <=>
    normalString delimiter;

    if
    (
        str.find('<') != std::string::npos
     && str.find('=') != std::string::npos
     && str.find('>') != std::string::npos
    )
    {
        delimiter = "<=>";

        //- forward and backward reaction is used
        //  prevent due to reactants and products
        if (i == 0)
        {
            kfkb_.push_back(1);
        }
    }
    else if
    (
        str.find('=') != std::string::npos
     && str.find('>') != std::string::npos
    )
    {
        delimiter = "=>";

        //- only forward reaction is used
        //  prevent due to reactants and products
        if (i == 0)
        {
            kfkb_.push_back(0);
        }
    }
    else if
    (
        str.find('=') != std::string::npos
    )
    {
        delimiter = "=";

        //- forward and backward reaction is used
        //  prevent due to reactants and products
        if (i == 0)
        {
            kfkb_.push_back(1);
        }
    }
    //- something went wrong
    else
    {
        std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
            << __LINE__ << " ++ " << str << std::endl;
        std::terminate();
    }

    //- return reactants
    if (i == 0)
    {
        return str.substr(0, str.find(delimiter));
    }
    //- return products
    else if (i == 1)
    {
        return str.substr(str.find(delimiter)+delimiter.size(), str.size());
    }
    //- i is not defined
    else
    {
        std::cerr<< " ++ ERROR in " << __FILE__ << " line no. "
            << __LINE__ << " ++ i is not defined. i = " << i << std::endl;
        std::terminate();
    }
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

    unsigned int kb{0};
    unsigned int fO{0};
    unsigned int lP{0};

    //- get backward reactions
    forAll(kfkb_, i)
    {
        if (kfkb_[i] == 1)
        {
            kb++;
        }
    }

    //- get fall off reactions
    forAll(nTROE_, i)
    {
        if (nTROE_[i] == 1)
        {
            fO++;
        }
    }

    //- get low pressure reactions
    forAll(nLOW_, i)
    {
        if (nLOW_[i] == 1)
        {
            lP++;
        }
    }

//    std::cout<< "----------------------------------------------------------\n"
//             << "                   CHEMISTRY SUMMARY                      \n"
//             << "----------------------------------------------------------\n"
//             << std::left
//             << std::setw(40) << " No. of elements: "
//             << std::setw(10) << elements_.size() << "\n"
//             << std::setw(40) << " No. of species in reactions: "
//             << std::setw(20) << species_.size() << "\n"
//             << std::setw(40) << " No. of elementar reactions: "
//             << std::setw(20) << n_ << "\n"
//             << std::setw(40) << " No. of low pressure reactions: "
//             << std::setw(20) << lP << "\n"
//             << std::setw(40) << " No. of fall off reactions: "
//             << std::setw(20) << fO << "\n"
//             << std::setw(40) << " No. of irreversible reactions: "
//             << std::setw(20) << kfkb_.size() - kb << "\n"
//             << std::setw(40) << " No. of reversible reactions: "
//             << std::setw(20) << kb << "\n"
//             << std::setw(40) << " No. of dublicated reactions: "
//             << std::setw(20) << nDuplicate_ << "\n"
//             << "----------------------------------------------------------\n"
//             << std::setw(40) << " Matrix size nu_: "
//             << std::setw(2) << nu_.size() << "x"
//             << nu_[0].size() << "\n"
//             << std::setw(40) << " Matrix size Mvalue_: "
//             << std::setw(2) << Mvalue_.size() << "x" << Mvalue_[0].size() << "\n"
//             << std::setw(40) << " Matrix size ArrheniusCoeffs_: "
//             << std::setw(2) << arrheniusCoeffs_.size() << "x"
//             << arrheniusCoeffs_[0].size() << "\n"
//             << std::setw(40) << " Matrix size TROECoeffs_: "
//             << std::setw(2) << TROECoeffs_.size() << "x"
//             << TROECoeffs_[0].size() << "\n"
//             << std::setw(40) << " Matrix size LOWCoeffs_: "
//             << std::setw(2) << LOWCoeffs_.size() << "x"
//             << LOWCoeffs_[0].size() << "\n"
//             << std::setw(40) << " Vector size kfkb_: "
//             << std::setw(2) << kfkb_.size() << "\n"
//             << "----------------------------------------------------------\n";
            std::cout << std::left << "reaction no." << std::setw(45) << " " << "TBR\tENHANCED\tLOW\tTROE\tSRI\n";
             for (unsigned int i=0; i<n_; i++)
             {
                std::cout << std::left << i+1 << ": " << std::setw(55) << elementarReaction_[i]
                << TBR_[i] << "\t"
                << ENHANCE_[i] << "\t"
                << LOW_[i] << "\t"
                << TROE_[i] << "\t"
                << SRI_[i] << "\t"
                << "\n";
             }
}
