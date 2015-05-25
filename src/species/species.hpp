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
#include <string>

//- user def. headers
#include "../definitions/typedef.hpp"
#include "../thermodynamic/thermodynamic.hpp"

#ifndef SPECIES_HPP
#define SPECIES_HPP


class Species : public Thermodynamic
{
    public:

        //- Constructor
        Species();
        Species(normalString);

        //- Destructor
        ~Species();


    public:

        //- Fuctions

            //- set the name of species
            void setName(const normalString&);

            //- set the molecular weight of species
            void setMW(const scalar);

            //- return name of species
            const normalString name();

            //- retrun MW [mol/kg]
            const scalar MW();


    private:

        //- Name
        normalString name_;

        //- Molecular weight [kg/mol]
        scalar MW_;

        //- Mass fraction [-]
        scalar Y_;

        //- Mol fraction [-]
        scalar X_;

        //- Concentration [-]
        scalar Con_;
};

#endif // SPECIES_HPP
