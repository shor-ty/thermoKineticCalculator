/*---------------------------------------------------------------------------*\
| =========                |                                                  |
| \\      /  A utomatic    | Holzmann-cfd                                     |
|  \\    /   F lamelet     | Version: 1.0                                     |
|   \\  /    C onstructor  | Web: www.Holzmann-cfd.de                         |
|    \\/                   |                                                  |
\*---------------------------------------------------------------------------*/
/*
»
» Description:
»   This class contains information about the elements:
»       + Name
»       + Molecular weight
»
»
» Used:
»   For calculating the net rate of Elements
»
\*---------------------------------------------------------------------------*/
//- system headers

//- user def. headers
#include "mixtureFraction.hpp"

//- constructor with species vector
MixtureFraction::MixtureFraction
(
    const stringField& species,
    const Chemistry* chemistryObj
)
{
    species_ = species;
    pChemistry = chemistryObj;

    //- increase all vectors
    forAll(species_, i)
    {
        //- mass fraction increment
        Y_.push_back(0);
        Yf_.push_back(0);
        Yo_.push_back(0);

        //- mole fraction increment
        X_.push_back(0);
        Xf_.push_back(0);
        Xo_.push_back(0);

        //- is fuel or oxidizer bool vector
        fuel_.push_back(false);
        oxidizer_.push_back(false);
    }
}


//- destructor
MixtureFraction::~MixtureFraction()
{
}

//- read afcDict
void MixtureFraction::readAFCDict
(
    const normalString& fileName
)
{
    //- output
    std::cout << "Reading AFCDict (" << fileName << ")\n\n";

    //- read the whole file and store it
    stringField fileContent = Transport::openFile(fileName);

    bool fuelFound{false};
    bool oxidizerFound{false};
    bool fuelOxidizerFound{false};
    bool inDictionary{false};

    //- loop through the fileContent
    forAll(fileContent, line)
    {
        //- STEP 1: remove all whitespaces
        stringField tmp = Transport::splitString(fileContent[line]);

        //- if not empty, and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            //- moleFractionFuel and moleFractionOxidizer
            if
            (
                tmp[0] == "moleFractionFuel"
             || tmp[0] == "moleFractionOxidizer"
            )
            {
                if (tmp[0] == "moleFractionFuel")
                {
                    fuelFound = true;
                }
                else if (tmp[0] == "moleFractionOxidizer")
                {
                    oxidizerFound = true;
                }

                fuelOxidizerFound = true;
            }
            //- fuel and oxidizer dictionary
            else if
            (
                fuelOxidizerFound
             && !inDictionary
            )
            {
                if (tmp[0] == "{")
                {
                    inDictionary = true;
                }
            }
            //- set fuel and oxidizer boundary conditions
            else if
            (
                fuelOxidizerFound
             && inDictionary
            )
            {
                if (tmp[0] == "}")
                {
                    inDictionary = false;
                    fuelOxidizerFound = false;
                    oxidizerFound = false;
                    fuelFound = false;
                }
                else
                {
                    //- oxidizer boundary
                    if (oxidizerFound)
                    {
                        bool speciesFound{false};

                        //- find species
                        forAll(species_, i)
                        {
                            if (species_[i] == tmp[0])
                            {
                                //- set oxidizer parameters
                                Xo_[i] = stod(tmp[1]);
                                oxidizer_[i] = true;
                                speciesFound = true;
                            }
                        }

                        //- check if species found
                        if (!speciesFound)
                        {
                            std::cerr<< "\n ++ ERROR: Species " << tmp[0]
                                     << " not found in chemistry class."
                                     << __FILE__ << " line " << __LINE__
                                     << std::endl;
                            std::terminate();
                        }
                    }
                    //- fuel boundary
                    else if (fuelFound)
                    {
                        bool speciesFound{false};

                        //- find species
                        forAll(species_, i)
                        {
                            if (species_[i] == tmp[0])
                            {
                                //- set fuel parameters
                                Xf_[i] = stod(tmp[1]);
                                fuel_[i] = true;
                                speciesFound = true;
                            }
                        }

                        //- check if species found
                        if (!speciesFound)
                        {
                            std::cerr<< "\n ++ ERROR: Species " << tmp[0]
                                     << " not found in chemistry class."
                                     << __FILE__ << " line " << __LINE__
                                     << std::endl;
                            std::terminate();
                        }
                    }
                }
            }

            //- set oxidizer temperature
            else if
            (
                tmp[0] == "temperatureOxidizer"
            )
            {
                To_ = stod(tmp[1]);
            }

            //- set fuel temperature
            else if
            (
                tmp[0] == "temperatureFuel"
            )
            {
                Tf_ = stod(tmp[1]);
            }
        }
    }
}


//- return species i
const stringField MixtureFraction::species() const
{
    return species_;
}


//- calculate mass fraction out of mole fraction
void MixtureFraction::XToY()
{
    forAll(species_, i)
    {
        Yf_ = pChemistry->MW(i) /

             species[id].MW() /
//                    MMW_fuel *
//                    species[id].Xf()
    }
}
