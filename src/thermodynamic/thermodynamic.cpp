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
#include "thermodynamic.hpp"

//- constructor
Thermodynamic::Thermodynamic()
:
    thermoBool(false)
{
}

//- destructor
Thermodynamic::~Thermodynamic()
{
}


void Thermodynamic::readThermodynamicFile
(
    const normalString& fileName,
    const stringField& species
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
                for (;line < fileContent.size();line++)
                {
                    tmp = splitString(fileContent[line]);

                    //- if line is not empty and no comment, proceed
                    if
                    (
                        !tmp.empty()
                     && tmp[0][0] != '!'
                    )
                    {
                        readNASA
                        (
                            fileContent,
                            line,
                            species
                        );
                    }
                }
            }
        }
    }
}


void Thermodynamic::readNASA
(
    const stringField& fileContent,
    const unsigned int& line,
    const stringField& species
)
{
    //- split line with delimiter ' '
    stringField tmp = splitString(fileContent[line]);


    //- LINE 1
    if (fileContent[line][79] == '1')
    {
        //- species name
        stringField speciesInThermo = splitString
            (
                fileContent[line].substr(0,18)
            );

        //- species formula
        normalString atomicComposition = fileContent[line].substr(24,19);

        //- get ID of species saved in chemistry class
        forAll(species, id)
        {
            if (species[id] == speciesInThermo[0])
            {
                NASA_[id] = true;
                speciesId_ = id;

                //- calculate molecular weight
                calcMolecularWeight
                (
                    atomicComposition
                );
            }
        }
    }
}


void Thermodynamic::thermodynamicDataIncrement()
{
    //- increment matrix of NASACoeffs
    NASACoeffs_.push_back(std::vector<double>(14));

    //- increment bool vector if NASACoeffs found for species
    NASA_.push_back(false);

    //- increment vector of molecular weight of species
    MW_.push_back(0);
}


void Thermodynamic::calcMolecularWeight
(
    const normalString& composition
)
{
    stringField tmp = splitString(composition);
    normalString composition_;

    //- put all together
    forAll(tmp, i)
    {
        composition_ += tmp[i];
    }

    bool foundLetter{false};
    bool foundNumber{false};

    bool lastNumber{false};

    normalString letter;
    normalString tmp2;

    scalar multiplicator;
    scalar tmpMW{0};


   //- loop through all single letter of the species
    for (unsigned int pos = 0; pos < composition_.size(); pos++)
    {
        foundLetter = false;
        foundNumber = false;

        //- compare single letter with ASCII table (LETTERS)
        for (unsigned int j = 65; j <= 90; j++)
        {
            char c = static_cast<char>(j);

            if
            (
                c == composition_[pos]
            )
            {
                foundLetter = true;
                foundNumber = false;
                letter = c;
                break;
            }
            else
            {
                foundNumber = true;
                foundLetter = false;
            }
        }

        //- if letter found but befor there was a number (new species)
        if (foundLetter && lastNumber)
        {
            //- calculate MW for atomic element
            tmpMW += atomicWeight
            (
                tmp2,
                multiplicator
            );

            //- reset the string
            tmp2.clear();
            multiplicator = 0;
            lastNumber = false;
        }

        //- if a letter is found
        if (foundLetter)
        {
            //- last letter
            tmp2 += composition_[pos];
        }


        if (foundNumber)
        {
            //- last number
            lastNumber = true;

            multiplicator += atof(composition_.substr(pos,1).c_str());
        }

        //- if pos is at last position
        if (pos == composition_.size()-1)
        {
            //- calculate MW for atomic element
            tmpMW += atomicWeight
            (
                tmp2,
                multiplicator
            );

            //- save
            MW_[speciesId_] = tmpMW;
        }
        //- todo | logfile of Molecular weight
    }
}


bool Thermodynamic::NASA
(
    const unsigned int& i
) const
{
    return NASA_[i];
}


scalar Thermodynamic::MW
(
    const unsigned int& i
) const
{
    return MW_[i];
}

//- set phase
void Thermodynamic::setPhase( const normalString& phase_ )
{
    phase = phase_;
}

//- set low temperature
void Thermodynamic::setPolyTemperature
(
    const scalar& lowTemp_,
    const scalar& comTemp_,
    const scalar& higTemp_
)
{
    lowTemp = lowTemp_;
    comTemp = comTemp_;
    higTemp = higTemp_;
}
