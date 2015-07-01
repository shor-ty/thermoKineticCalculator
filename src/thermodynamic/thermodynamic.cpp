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

    bool found{false};
    bool lastLetter{false};
    normalString letter;
    normalString multiplicator;
    normalString tmp2;

   //- loop through all single letter of the species
    for (unsigned int pos = 0; pos < composition_.size(); pos++)
    {
        found = false;
        //- compare single letter with ASCII table (LETTERS)
        for (unsigned int j = 65; j <= 90; j++)
        {
            char c = static_cast<char>(j);

            if
            (
                c == composition_[pos]
            )
            {
                found = true;
                lastLetter = true;
                letter = c;
                break;
            }
        }

        if
        (
            found
         && lastLetter
        )
        {
            tmp2 += letter;
        }


        if (!found)
        {
            lastLetter = false;
           multiplicator += composition_[pos];
        }
        std::cout << composition << ">>>> " << tmp2 << " -> " << multiplicator << "\n";

        if (!found) tmp2.clear();
    }
}

stringField Thermodynamic::openFile
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


stringField Thermodynamic::splitString
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


stringField Thermodynamic::splitString
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


//- functions
//-----------

    //- set NASA coefficents
//    void Thermodynamic::setPolyCoeffs(const scalarField& polyCoeffs_)
//    {
//        polyCoeffs = polyCoeffs_;
//    }
//
//    //- calculate heat capacity cp [J/(molK)]
//    scalar Thermodynamic::cp(const scalar T_) const
//    {
//        int i{0};
//
//        //- get information of temperature range
//        whichPolyCoeffs(T_, i);
//
//        // return heat capacity [J/(molK)]
//        return ((   polyCoeffs[i]
//                  + polyCoeffs[i+1]*T_
//                  + polyCoeffs[i+2]*pow(T_,2)
//                  + polyCoeffs[i+3]*pow(T_,3)
//                  + polyCoeffs[i+4]*pow(T_,4) )
//                  * R);
//    }
//
//    //- calculate enthalpy [J/mol]
//    scalar Thermodynamic::h(const scalar T_) const
//    {
//        int i{0};
//
//        //- get information of temperature range
//        whichPolyCoeffs(T_, i);
//
//        // return enthalpy [J/(mol)]
//        return ((   polyCoeffs[i]
//                  + polyCoeffs[i+1]*T_/2
//                  + polyCoeffs[i+2]*pow(T_,2)/3
//                  + polyCoeffs[i+3]*pow(T_,3)/4
//                  + polyCoeffs[i+4]*pow(T_,4)
//                  + polyCoeffs[i+5]/T_)
//                  *T_*R);
//    }
//
//    //- return ref enthalpy h0 [J/mol]
//    scalar Thermodynamic::h0() const
//    {
//        int i{0};
//
//        //- get information of temperature range
//        whichPolyCoeffs(298, i);
//
//        // return ref enthalpy h0 [J/mol]
//        return polyCoeffs[i+5]*R;
//    }
//
//    //- calculate entropy [J/(molK)]
//    scalar Thermodynamic::s(const scalar T_) const
//    {
//        int i{0};
//
//        //- get information of temperature range
//        whichPolyCoeffs(T_, i);
//
//        // return entropy [J/(molK)]
//        return ((   polyCoeffs[i]*log(T_)
//                  + polyCoeffs[i+1]*T_
//                  + polyCoeffs[i+2]*pow(T_,2)/2
//                  + polyCoeffs[i+3]*pow(T_,3)/3
//                  + polyCoeffs[i+4]*pow(T_,4)/4
//                  + polyCoeffs[i+6])
//                  *R);
//    }
//
//    //- set thermodynamic bool to true
//    void Thermodynamic::thermodynamicTrue()
//    {
//        thermoBool = true;
//    }
//
//    //- return thermodynamic status
//    bool Thermodynamic::thermodynamicStatus() const
//    {
//        return thermoBool;
//    }
//
//    //- check temperature range
//    void Thermodynamic::whichPolyCoeffs(const scalar& T_, int& i) const
//    {
//        //- FIXME EXTEND TO WHICH SPECIES
//        if(T_ < lowTemp)
//        {
//            std::cerr << "  ++ Warning: NASA-Polynomial is not defined for temperature  T=" << T_ << " K\n"
//                      << "  ++ Warning: Using low temperature coefficents for calculation...\n";
//            i = 7;
//        }
//        else if(T_ >= lowTemp && T_  <= comTemp)
//        {
//            i = 7;
//        }
//        else if(T_ > comTemp && T_  <= higTemp)
//        {
//            i = 0;
//        }
//        else
//        {
//            std::cerr << "  ++ Warning: NASA-Polynomial is not defined for temperature  T=" << T_ << " K\n"
//                     << "  ++ Warning: Using high temperature coefficents for calculation...\n";
//            i = 0;
//        }
//    }
