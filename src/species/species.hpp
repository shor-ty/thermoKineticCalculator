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
#include "../chemistry/chemistry.hpp"

#ifndef SPECIES_HPP
#define SPECIES_HPP


class Species
:
    public Thermodynamic,
    public Chemistry
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

            //- mol fraction

                //- set mol fraction X 0 < Z < 1
                void setX(const scalar&);

                //- set mol fraction X in pure fuel
                void setXf(const scalar&);

                //- set mol fraction X in pure oxidizer
                void setXo(const scalar&);

                //- return mol fraction X 0 < Z < 1
                const scalar X() const;

                //- return mol fraction X in pure fuel
                const scalar Xf() const;

                //- return mol fraction X in pure oxidizer
                const scalar Xo() const;


            //- mass fraction

                //- set mass fraction Y 0 < Z < 1
                void setY(const scalar&);

                //- set mass fraction Y in pure fuel
                void setYf(const scalar&);

                //- set mass fraction Y in pure oxidizer
                void setYo(const scalar&);

                //- return mass fraction Y 0 < Z < 1
                const scalar Y() const;

                //- return mass fraction Y in pure fuel
                const scalar Yf() const;

                //- return mass fraction Y in pure oxider
                const scalar Yo() const;


            //- temperature

                //- set temperature of pure fuel
                void setTf(const scalar&);

                //- set temperature of pure oxidizer
                void setTo(const scalar&);

                //- get temperature of pure fuel
                const scalar Tf() const;

                //- get temperature of pure oxidizer
                const scalar To() const;


            //- return name of species
            const normalString name() const;

            //- retrun MW [g/mol]
            const scalar MW() const;



    private:

        //- name
        normalString name_;

        //- molecular weight [g/mol]
        scalar MW_;


        //- mass fraction

            //- mass fraction [-]
            scalar Y_;

            //- mass fraction in pure fuel [-]
            scalar Yf_;

            //- mass fraction in pure oxidizer [-]
            scalar Yo_;


        //- mol fraction

            //- mol fraction [-]
            scalar X_;

            //- mol fraction in pure fuel [-]
            scalar Xf_;

            //- mol fraction in pure oxidizer [-]
            scalar Xo_;


        //- temperature [K]
        scalar T_;

        //- temperature of pure fuel
        scalar Tf_;

        //- temperature of pure oxidizer
        scalar To_;

        //- concentration [-]
        scalar Con_;

        //- is part of fuel mixture
        bool fuel_;

        //- is part of oxidizer mixture
        bool oxidizer_;
};

#endif // SPECIES_HPP
