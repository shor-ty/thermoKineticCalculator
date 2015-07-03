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
#include "transport.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Transport::Transport()
{
}

Transport::~Transport()
{
}


void Transport::transportDataIncrement()
{
    //- increment all vectors that contain transport stuff

        //- vector of Leonard-Jones Potential
        LJP_.push_back(0);

        //- vector of Leonard-Jones collision diameter
        LJCD_.push_back(0);

        //- vector of dipole moment
        DPM_.push_back(0);

        //- vector of polarizability
        P_.push_back(0);

        //- vector of rotational relaxation collision number
        RRCN_.push_back(0);

        //- vector for transport properties (found or not)
        TRANS_.push_back(false);
}


void Transport::readTransportFile
(
    const normalString& fileName,
    const stringField& species
)
{
    //- output
    std::cout << "Reading transport properties (" << fileName << ")\n";

    stringField fileContent = openFile(fileName);

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
            //  find the species and save the transport data
            forAll(species, id)
            {
                if (species[id] == tmp[0])
                {
                    //- bool to check if all species got transport data
                    TRANS_[id] = true;

                    //- transport data
                    LJP_[id] = stoi(tmp[1]);
                    LJCD_[id] = stod(tmp[2]);
                    DPM_[id] = stod(tmp[3]);
                    P_[id] = stod(tmp[4]);
                    RRCN_[id] = stod(tmp[5]);
                }
            }
        }
    }
}


stringField Transport::openFile
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


stringField Transport::splitString
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


stringField Transport::splitString
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


const bool Transport::TRANS
(
    const unsigned int& i
) const
{
    return TRANS_[i];
}
