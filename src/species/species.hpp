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

            //- set fuel | set oxidizer
            void setFuel();
            void setOxidizer();

            //- return fuel | oxidizer
            const bool fuel() const;
            const bool oxidizer() const;

            //- set mol fraction X
            void setX(const scalar&);

            //- return mol fraction X
            const scalar X() const;

            //- return name of species
            const normalString name() const;

            //- retrun MW [mol/kg]
            const scalar MW() const;



    private:

        //- name
        normalString name_;

        //- molecular weight [kg/mol]
        scalar MW_;

        //- mass fraction [-]
        scalar Y_;

        //- mol fraction [-]
        scalar X_;

        //- concentration [-]
        scalar Con_;

        //- is part of fuel mixture
        bool fuel_;

        //- is part of oxidizer mixture
        bool oxidizer_;
};

#endif // SPECIES_HPP
