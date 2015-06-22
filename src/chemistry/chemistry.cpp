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

//- user def. headers
#include "chemistry.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Chemistry::Chemistry()
:
    r_(0),
    rDuplicate_(0)
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
                    rDuplicate_++;

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
                        //- re-arange the line
                        //  needed if reactions are not together e.g.
                        //  H2+     O2  =  OH+OH
                        //  need that form
                        //  H2+O2=OH+OH
                        normalString reaction;

                        //- add all strings exept the last three
                        //  last three are arrhenius coefficients
                        for (unsigned int i=0; i<lineArray.size()-3; i++)
                        {
                            reaction += lineArray[i];
                        }

                        //- temp vector which contains the arrhenius coefficient
                        scalarField arrhenius;

                        //- arrhenius coefficients
                        for (unsigned int i=lineArray.size()-3; i < lineArray.size(); i++)
                        {
                            arrhenius.push_back(stod(lineArray[i]));
                        }

                        //- insert the reaction to the vector
                        //  added as string
                        elementarReaction_.push_back(reaction);

                        //- increment matrix nu
                        nu_.push_back(std::vector<double>(species_.size()));

                        //- increment matrix M
                        M_.push_back(std::vector<double>(species_.size()));

                        //- increment arrhenius, TROE and LOW matrix
                        arrheniusCoeffs_.push_back(std::vector<double>(3));
                        TROECoeffs_.push_back(std::vector<double>(3));
                        LOWCoeffs_.push_back(std::vector<double>(3));

                        //- increment fall off vector
                        fO_.push_back(0);

                        //- increment low pressure vector
                        lP_.push_back(0);

                        //- splitting the reaction into reactants and products
                        //  additionall the matrix kfkb is modified
                        //  -  only forward reaction kfkb := 0
                        //  -  both reactions are used    := 1
                        //  DONT MESS IT UP WITH THE FOLLOWING FUNCTION:
                        //  meaning:
                        //  -  splitReaction(string, 0): return reactants
                        //  -  splitReaction(string, 1): return products
                        normalString reactants = splitReaction(reaction, 0);
                        normalString products = splitReaction(reaction, 1);

                        //- use the reaction and get the following points
                        //  get stochiometric factors for each species
                        //  set if LOW and TROE
                        //  set arrhenius coefficients
                        //  set LOW and TROE coefficients
                        //  set all other stuff
                        //  function:
                        //  +   reactants | products
                        //  +   0 | 1 for switching between reactants and products
                        //  +   coefficients for that reaction
                        //  +   fileContent and line number for LOW and TROE
                        updateAllMatrix(reactants, 0, arrhenius, fileContent, line);
                        updateAllMatrix(products, 1, arrhenius, fileContent, line);


                        //- increment reaction counter
                        r_++;
                    }
                }
            }
        //- REACTIONS end
        }
    //- loop through the file content end
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


stringField Chemistry::splitString
(
    const normalString str
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
    const normalString& str,
    const unsigned int i
)
{
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
        if (kfkb_.size() != r_)
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
        if (kfkb_.size() != r_)
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
        if (kfkb_.size() != r_)
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


void Chemistry::updateAllMatrix
(
    const normalString& reaction,
    const unsigned int i,
    const scalarField& arrhenius,
    const stringField& fileContent,
    const unsigned int& line
)
{
    //- EXCEPTION 1
    //  check (+M)
    normalString pattern = "(+M)";
    std::size_t found = reaction.find(pattern);

    //- THIRD BODY REACTION WITH (+M)
    //  LOW and TROE
    if (found != std::string::npos)
    {
        //- modify the reaction (remove (+M))
        normalString modifiedReaction = reaction.substr(0,found);

        //- spilt into the single species using delimiter '+'
        std::stringstream tmp(modifiedReaction);
        normalString element;
        stringField elements;

        while (std::getline(tmp, element, '+'))
        {
            elements.push_back(element);
        }


        //- search the ID of species
        forAll(species_, id)
        {
            forAll(elements, elem)
            {
                if (species_[id] == elements[elem])
                {
                    int nuTmp_ = nu_[r_][id];

                    //- sign of nu depend on the side
                    //  Definition
                    //  +  reactants := negativ
                    //  +  products := positiv
                    if (i == 0)
                    {
                        nuTmp_--;
                    }
                    else if (i == 1)
                    {
                        nuTmp_++;
                    }
                    //- something went wrong
                    else
                    {
                        std::cerr<< " ++ ERROR in " << __FILE__
                            << " line no." << __LINE__ << " ++ i is "
                            << "not defined. i = " << i << ". Species"
                            << " not defined: " << elements[elem]
                            << std::endl;
                        std::terminate();
                    }

                    //- set the stochiometric factors
                    nu_[r_][id] = nuTmp_;


                    //- set arrhenius coefficients
                    //  +   0 := A0
                    //  +   1 := b
                    //  +   2 := Ea
                    arrheniusCoeffs_[r_][0] = arrhenius[0];
                    arrheniusCoeffs_[r_][1] = arrhenius[1];
                    arrheniusCoeffs_[r_][2] = arrhenius[2];

                    //- check if LOW
                    normalString line1 = fileContent[line+1];
                    stringField  line1Array = splitString(line1);

                    if (line1Array[0] == "LOW/")
                    {
                        //- set low pressure (LOW)
                        lP_[r_] = 1;

                        //- next check
                        //  TROE used?
                        normalString line2 = fileContent[line+2];
                        stringField  line2Array = splitString(line2);


                        if (line2Array[0] == "TROE/")
                        {
                            //- set fall off (TROE)
                            fO_[r_] = 1;

                            //- set TROES formula coefficient
                            TROECoeffs_[r_][0] = stod(line2Array[1]);
                            TROECoeffs_[r_][1] = stod(line2Array[2]);
                            TROECoeffs_[r_][2] = stod(line2Array[3]);

                            //- set LOW pressure coefficient
                            LOWCoeffs_[r_][0] = stod(line1Array[1]);
                            LOWCoeffs_[r_][1] = stod(line1Array[2]);
                            LOWCoeffs_[r_][2] = stod(line1Array[3]);

                            //- check next
                            //  modified M?
                            normalString line3 = fileContent[line+3];

                            //- if modified M, next line does not contain '='
                            //  and contain '/'
                            if
                            (
                                line3.find('=') == std::string::npos
                             && line3.find('/') != std::string::npos
                            )
                            {
                                //- split the line with delimiter '/'
                                std::stringstream tmp(line3);
                                normalString element;
                                stringField elements;

                                while (std::getline(tmp, element, '/'))
                                {
                                    elements.push_back(element);
                                }

                                forAll(elements, elem)
                                {
                                    forAll(species_, id)
                                    {
                                        if (species_[id] == elements[elem])
                                        {
                                            M_[r_][id] =
                                                stod(elements[elem+1]);
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            //- set low pressure (LOW)
                            lP_[r_] = 1;

                            //- set LOW pressure coefficient
                            LOWCoeffs_[r_][0] = stod(line1Array[1]);
                            LOWCoeffs_[r_][1] = stod(line1Array[2]);
                            LOWCoeffs_[r_][2] = stod(line1Array[3]);

                            //- check next
                            //  modified M?
                            normalString line3 = fileContent[line+3];

                            //- if modified M, next line does not contain '='
                            //  and contain '/'
                            if
                            (
                                line3.find('=') == std::string::npos
                             && line3.find('/') != std::string::npos
                            )
                            {
                                //- split the line with delimiter '/'
                                std::stringstream tmp(line3);
                                normalString element;
                                stringField elements;

                                while (std::getline(tmp, element, '/'))
                                {
                                    elements.push_back(element);
                                }

                                forAll(elements, elem)
                                {
                                    forAll(species_, id)
                                    {
                                        if (species_[id] == elements[elem])
                                        {
                                            M_[r_][id] =
                                                stod(elements[elem+1]);
                                        }
                                    }
                                }
                            //- LOW - M modification end
                            }
                        //- LOW end
                        }
                    //- LOW TROE end
                    }
                //- if condition if species in reaction == species
                }
            //- loop over species which are in reaction
            }
        //- loop over all species which involved in elementar reaction
        }
    //- THIRD BODY REACTION TROE end (EXCEPTION 1)
    }
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
    forAll(fO_, i)
    {
        if (fO_[i] == 1)
        {
            fO++;
        }
    }

    //- get low pressure reactions
    forAll(lP_, i)
    {
        if (lP_[i] == 1)
        {
            lP++;
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
             << std::setw(20) << r_ << "\n"
             << std::setw(40) << " No. of low pressure reactions: "
             << std::setw(20) << lP << "\n"
             << std::setw(40) << " No. of fall off reactions: "
             << std::setw(20) << fO << "\n"
             << std::setw(40) << " No. of forward reactions: "
             << std::setw(20) << kfkb_.size() << "\n"
             << std::setw(40) << " No. of backward reactions: "
             << std::setw(20) << kb << "\n"
             << std::setw(40) << " No. of dublicated reactions: "
             << std::setw(20) << rDuplicate_ << "\n"
             << "----------------------------------------------------------\n"
             << std::setw(40) << " Matrix size nu_: "
             << std::setw(2) << nu_.size() << "x"
             << nu_[0].size() << "\n"
             << std::setw(40) << " Matrix size M_: "
             << std::setw(2) << M_.size() << "x" << M_[0].size() << "\n"
             << std::setw(40) << " Matrix size ArrheniusCoeffs_: "
             << std::setw(2) << arrheniusCoeffs_.size() << "x"
             << arrheniusCoeffs_[0].size() << "\n"
             << std::setw(40) << " Matrix size TROECoeffs_: "
             << std::setw(2) << TROECoeffs_.size() << "x"
             << TROECoeffs_[0].size() << "\n"
             << std::setw(40) << " Matrix size LOWCoeffs_: "
             << std::setw(2) << LOWCoeffs_.size() << "x"
             << LOWCoeffs_[0].size() << "\n"
             << std::setw(40) << " Vector size kfkb_: "
             << std::setw(2) << kfkb_.size() << "\n"
             << "----------------------------------------------------------\n";

             forAll(kfkb_, i)
             std::cout << i << ": " << kfkb_[i] << " --> " << elementarReaction_[i] << "\n";
}
