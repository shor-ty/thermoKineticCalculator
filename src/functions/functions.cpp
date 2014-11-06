/*---------------------------------------------------------------------------*\
| =========                |                                                  |
| \\      /  A utomatic    | Holzmann-cfd                                     |
|  \\    /   F lamelet     | Version: 1.0                                     |
|   \\  /    C constructor | Web: www.Holzmann-cfd.de                         |
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
#include <iostream>
#include <regex>
#include <fstream>
#include <string>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

//- user def. headers
#include "functions.hpp"
#include "species.hpp"
#include "chemistry.hpp"

//- function definition
//---------------------

    //- check input parameters
    void checkInputParameters( const int & noOfArguments, const char* argument1, const char* argument2)
    {
        if( noOfArguments > 7 )
        {
            showErrorArgumentInput();
            std::terminate();
        }
        else if( noOfArguments < 7)
        {
            showErrorArgumentInput();
            std::terminate();
        }
        else
        {
            std::string arg1 = argument1;
            std::string arg2 = argument2;

            std::regex kin ("(.*)(\\.kin)");
            std::regex tdc ("(.*)(\\.tdc)");
            std::regex tra ("(.*)(\\.tra)");

            if( !(( arg1 == "-kinetic" && std::regex_match( arg2,kin )) || ( arg1 == "-thermo" && std::regex_match( arg2,tdc )) || ( arg1 == "-transport" && std::regex_match( arg2,tra ))) )
            {
                showErrorWrongInput();
                std::terminate();
            }
        }
    }

    //- show start info
    void showInfo()
    {
        std::cout <<    "/*---------------------------------------------------------------------------*\\\n" <<
                        "| =========                |                                                  |\n" <<
                        "| \\\\      /  A utomatic    | Holzmann-cfd                                     |\n" <<
                        "|  \\\\    /   F lamelet     | Version: 1.0                                     |\n" <<
                        "|   \\\\  /    C constructor | Web: www.Holzmann-cfd.de                         |\n" <<
                        "|    \\\\/                   |                                                  |\n" <<
                        "\\*---------------------------------------------------------------------------*/\n" << std::endl;
    }

    //- show error; wrong parameters
    void showErrorArgumentInput()
    {
        std::cerr << "\n ++ ERROR: Wrong amount of arguments; expected: ./flameletGeneration -kinetic *.kin -thermodynamic *.tdc -transport *.tra";
        std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__ << std::endl;
    }

    //- show error; wrong input
    void showErrorWrongInput()
    {
        std::cerr << "\n ++ ERROR: Wrong file endings or wrong arguments; expected: ./flameletGeneration -kinetic *.kin -thermodynamic *.tdc -transport *.tra";
        std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__ << std::endl;
    }

    //- create array of all species objects
    std::vector<Species> createSpeciesObjects( const std::string& file_path )
    {
        std::ifstream file_tra;

        //-
        file_tra.open( file_path.c_str(), std::ios::in );
        if( !file_tra.good() )
        {
            std::cerr << "\n ++ ERROR: Could not open file \"" << file_path << "\" abort...";
            std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__ << std::endl;
            std::terminate();
        }

        //- temp variables
        std::string fileLine;
        int nSpecies{0};

        file_tra.seekg( 0L, std::ios::beg );

        //- how many species
        while ( !file_tra.eof() )
        {
            std::getline( file_tra, fileLine );
            nSpecies++;
        }
        file_tra.close();

        //- return dynamic object
        return std::vector<Species>(nSpecies-1);
    }

    //- get transport data
    void getTransportData( const std::string& file_path, std::vector<Species>& vecOfAllSpecies )
    {
        std::ifstream file_tra;

        //-
        file_tra.open( file_path.c_str(), std::ios::in );
        if( !file_tra.good() )
        {
            std::cerr << "\n ++ ERROR: Could not open file \"" << file_path << "\" abort...";
            std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__ << std::endl;
            std::terminate();
        }

        //- actual readed line
        std::string fileLine;
        //- actual line number
        int i{0};
        bool lessArg{false};

        //- read line per line and insert values to objects
        while ( !file_tra.eof() )
        {
            std::getline( file_tra, fileLine );
            boost::char_separator<char> sep(" ");
            boost::tokenizer<boost::char_separator<char>> elementArray( fileLine, sep );
            std::string::size_type sz;

            int j{1};
            for( const auto& argu : elementArray )
            {
                //- Check if too less arguments in line befor
                if ( lessArg && j == 1 )
                {
                    std::cerr << "\n ++ ERROR: transport.tra file is not in a correct form (too less arguments)";
                    std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__;
                    std::cerr << "\n ++ Line " << i << " is wrong in file '" << file_path << "'" << std::endl;
                    std::terminate();
                }
                //- First entry is name
                if ( j == 1 )
                {
                    vecOfAllSpecies[i].setID( i );
                    vecOfAllSpecies[i].setName( argu );
                    lessArg = true;
                }
                //- Second entry is geomtry
                else if ( j == 2 )
                {
                    vecOfAllSpecies[i].setGeometry( std::stoi( argu, &sz ) );
                    lessArg = true;
                }
                //- Third entry is potential
                else if ( j == 3 )
                {
                    vecOfAllSpecies[i].setLJPotential( std::stof( argu ) );
                    lessArg = true;
                }
                //- Fourth entry is collision diameter
                else if ( j == 4 )
                {
                    vecOfAllSpecies[i].setLJCollisionDiameter( std::stof( argu ) );
                    lessArg = true;
                }
                //- Fifth entry is dipol moment
                else if ( j == 5 )
                {
                    vecOfAllSpecies[i].setDipoleMoment( std::stof( argu ) );
                    lessArg = true;
                }
                //- Sixth entry is polarization
                else if ( j == 6 )
                {
                    vecOfAllSpecies[i].setPolarizability( std::stof( argu ) );
                    lessArg = true;
                }
                //- Fifth entry is rotation relaxation collision number
                else if ( j == 7 )
                {
                    vecOfAllSpecies[i].setRotRelaxCollNo( std::stof( argu ) );
                    lessArg = false;
                }
                else
                {
                    std::cerr << "\n ++ ERROR: transport.tra file is not in a correct form (too many arguments)";
                    std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__;
                    std::cerr << "\n ++ Line " << i+1 << " is wrong in file '" << file_path << "'" << std::endl;
                    std::terminate();
                }
                j++;
            }
            i++;
        }
        file_tra.close();
    }

    void getThermoData( const std::string& file_path, std::vector<Species>& vecOfAllSpecies )
    {
        std::ifstream file_tdc;

        //-
        file_tdc.open( file_path.c_str(), std::ios::in );
        if( !file_tdc.good() )
        {
            std::cerr << "\n ++ ERROR: Could not open file \"" << file_path << "\" abort...";
            std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__ << std::endl;
            std::terminate();
        }

        int l{0}, speciesID{-1}, lineNumber{1};
        //- read line per line and store data in objects
        while( !file_tdc.eof() )
        {
            std::string line;
            std::vector <std::string> elementInLine;

            //-
            std::getline( file_tdc, line);

            //- check if line is empty or not
            if ( !line.empty() )
            {
                split( elementInLine, line, boost::is_any_of( " " ), boost::token_compress_on );

                //- Check if first word in first line is "THERMO"
                if ( lineNumber == 0 && (elementInLine[0].compare( "THERMO" ) != 0) )
                {
                    std::cerr << "\n ++ ERROR: File \"" << file_path << "\" has wrong syntax in first line; abort...";
                    std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__;
                    std::cerr << "\n ++ Expected first word \"THERMO\" in first line " << lineNumber << std::endl;
                    std::terminate();
                }

                //- Start reading thermodynamic lines
                if (lineNumber >= 3)
                {

                    //- Every species has 4 lines
                    //- First line
                    if (l == 0)
                    {
                        //- Check if last entry takes number 1
                        if( elementInLine[elementInLine.size()-1].compare( "1" ) != 0 )
                        {
                            std::cerr << "\n ++ ERROR: '" << file_path << "' has wrong entrys; abort...";
                            std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__;
                            std::cerr << "\n ++ Expected '1' at the end of line " << lineNumber << " in file '" << file_path << "'" << std::endl;
                            std::terminate();
                        }

                        //- Search corresponding object
                        for( unsigned int i = 0; i < vecOfAllSpecies.size(); i++)
                        {
                            //- if species names are the same, set ID
                            if( vecOfAllSpecies[i].getName().compare( elementInLine[0] ) == 0 )
                            {
                                vecOfAllSpecies[i].setTransportThermodynamicToTrue();
                                speciesID = i;
                            }
                        }

                        //- just use species in thermodynamic
                        if( speciesID >= 0 )
                        {
                            //- Set low temperature
                            vecOfAllSpecies[speciesID].setLowTemperature( std::stoi( elementInLine[elementInLine.size()-4] ) );

                            //- Set mid temperature
                            vecOfAllSpecies[speciesID].setMidTemperature( std::stoi( elementInLine[elementInLine.size()-2] ) );

                            //- Set high temperatur
                            vecOfAllSpecies[speciesID].setHighTemperature( std::stoi( elementInLine[elementInLine.size()-3] ) );
                        }
                        //- l++ for line 2
                        l++;
                    }
                    //- FIRST LINE OF COEFFICENTS
                    else if ( l == 1 )
                    {
                        //- Check if last entry takes number 2
                        if( elementInLine[elementInLine.size()-1].compare( "2" ) != 0 )
                        {
                            std::cerr << "\n ++ ERROR: '" << file_path << "' has wrong entrys; abort...";
                            std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__;
                            std::cerr << "\n ++ Expected '2' at the end of line " << lineNumber << " in file '" << file_path << "'" << std::endl;
                            std::terminate();
                        }

                        //- just use species in thermodynamic
                        if( speciesID >= 0 )
                        {
                            //- Set the first line of coefficients
                            vecOfAllSpecies[speciesID].setFirstLineOfCoefficients
                            (
                                std::stod( line.substr( 0,15 ) ),
                                std::stod( line.substr( 15,15 ) ),
                                std::stod( line.substr( 30,15 ) ),
                                std::stod( line.substr( 45,15 ) ),
                                std::stod( line.substr( 60,15 ) )
                            );
                        }
                        else
                        {
                            std::cerr << "\n ** Warning: species in line " << lineNumber - 1<< " in file '" << file_path << "' is not in transportProperties...";
                        }


                        l++;
                    }
                    //- SECOND LINE OF COEFFICENTS
                    else if ( l == 2 )
                    {
                        if( elementInLine[elementInLine.size()-1].compare( "3" ) != 0 )
                        {
                            std::cerr << "\n ++ ERROR: '" << file_path << "' has wrong entrys; abort...";
                            std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__;
                            std::cerr << "\n ++ Expected '3' at the end of line " << lineNumber << " in file '" << file_path << "'" << std::endl;
                            std::terminate();
                        }

                        //- just use species in thermodynamic
                        if( speciesID >= 0 )
                        {
                            //- Set the second line of coefficients
                            vecOfAllSpecies[speciesID].setSecondLineOfCoefficients
                            (
                                std::stod( line.substr(0,15) ),
                                std::stod( line.substr(15,15) ),
                                std::stod( line.substr(30,15) ),
                                std::stod( line.substr(45,15) ),
                                std::stod( line.substr(60,15) )
                            );
                        }
                        l++;
                    }
                    //- THIRD LINE OF COEFFICENTS
                    else if ( l == 3 )
                    {
                        if( elementInLine[elementInLine.size()-1].compare( "4" ) != 0 )
                        {
                            std::cerr << "\n ++ ERROR: '" << file_path << "' has wrong entrys; abort...";
                            std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__;
                            std::cerr << "\n ++ Expected '4' at the end of line " << lineNumber << " in file '" << file_path << "'" << std::endl;
                            std::terminate();
                        }

                        //- just use species in thermodynamic
                        if( speciesID >= 0 )
                        {
                            //- Set the third line of coefficients
                            vecOfAllSpecies[speciesID].setThirdLineOfCoefficients
                            (
                                stod( line.substr(0,15) ),
                                stod( line.substr(15,15) ),
                                stod( line.substr(30,15) ),
                                stod( line.substr(45,15) )
                            );
                        }
                        l=0;
                        speciesID = -1;
                    }
                    else
                    {
                        std::cerr << "\n ++ ERROR: program problems; abort...";
                        std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__;
                        std::cerr << "\n ++ Variable 'l' has not defined values l valid  [0:3]" << std::endl;
                        std::terminate();
                    }
                }
            }
            //- rise line number
            lineNumber++;
        }
        file_tdc.close();
    }

    //- get kinetic data
    Chemistry generateAndGetKineticData( const std::string& file_path )
    {
        std::ifstream file_kin;

        //- generate chemistry object
        Chemistry chemObj;

        //-
        file_kin.open( file_path.c_str(), std::ios::in );
        if( !file_kin.good() )
        {
            std::cerr << "\n ++ ERROR: Could not open file \"" << file_path << "\" abort...";
            std::cerr << "\n ++ Error occur in file " << __FILE__ << " line " << __LINE__ << std::endl;
            std::terminate();
        }

        //- start count of elements, species and reactions
        bool startCountElements{false}, startCountSpecies{false}, startCountReactions{false};
        int countLinesElements{-1}, countLinesSpecies{-1}, countLinesReactions{-1};

        //- line number in file
        int lineNumber{1};

        //- read line per line and store data in objects
        while( !file_kin.eof() )
        {
            std::string line;
            std::vector <std::string> elementInLine;

            //-
            std::getline( file_kin, line);

            //- split line
            split( elementInLine, line, boost::is_any_of( " " ), boost::token_compress_on );

            //- element block
            {
                if( elementInLine[0].compare( "ELEMENTS" ) == 0 && !(startCountElements) )
                {
                    startCountElements = true;
                    countLinesElements = 0;
                }



                //- start counting elements after ELEMENTS till END
                if( startCountElements && elementInLine[0].compare( "END" ) != 0 )
                {
                    if ( countLinesElements > 0 )
                    {
                        chemObj.setElement( line );
                        countLinesElements++;
                    }
                    else
                    {
                        countLinesElements++;
                    }
                }
                //- if END reached;
                else
                {
                    //- reset variable
                    startCountElements = false;
                    countLinesElements = -1;
                }
            }

            //- species block
            {
                if( elementInLine[0].compare( "SPECIES" ) == 0 && !(startCountSpecies) )
                {
                    startCountSpecies = true;
                }

                //- start counting elements
                if( startCountSpecies && elementInLine[0].compare( "END" ) != 0 )
                {
                    //-
                    if(countLinesSpecies > 0)
                    {
                        chemObj.setSpecies( elementInLine );
                        countLinesSpecies++;
                    }
                    else
                    {
                        countLinesSpecies++;
                    }
                }
                //- if END reached;
                else
                {
                    //- reset variable
                    startCountSpecies = false;
                    countLinesSpecies = -1;
                }
            }
            //- reactions block
            {
                if( elementInLine[0].compare( "REACTIONS" ) == 0 && !(startCountReactions) )
                {
                    startCountReactions = true;
                }

                //- start counting elements
                if( startCountReactions && elementInLine[0].compare( "END" ) != 0 )
                {
                    //-
                    if(countLinesReactions > 0)
                    {
                        chemObj.setReactions( elementInLine );
                        countLinesReactions++;
                    }
                    else
                    {
                        countLinesReactions++;
                    }
                }
                //- if END reached;
                else
                {
                    //- reset variable
                    startCountReactions = false;
                    countLinesReactions = -1;
                }
            }

            lineNumber++;
        }
        file_kin.close();

        //- return object
        return chemObj;
    }
