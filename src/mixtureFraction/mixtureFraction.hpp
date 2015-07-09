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

#ifndef MixtureFraction_HPP
#define MixtureFraction_HPP


class MixtureFraction
{
    public:

        //- constructor

            //- constructor with argument of pointer of Chemistry obj
            //  and species vector
            MixtureFraction
            (
                const stringField&,
                const Chemistry*
            );


        //- destructor
        ~MixtureFraction();


    public:

        //- Pointer to chemistry class


        //- Fuctions

            //- read afcDict
            void readAFCDict
            (
                const normalString&
            );

            //- return species
            const stringField species() const;

            //- calculate mass fraction out of mole fraction
            void XToY();

            //- set fuel | set oxidizer
            void setFuel();
            void setOxidizer();

            //- return fuel | oxidizer
            const bool fuel() const;
            const bool oxidizer() const;

            //- show pointer adress
            void pointerAdress() const;

            //- mol fraction


            //- mass fraction


            //- temperature

                //- get temperature of pure fuel
                const scalar Tf() const;

                //- get temperature of pure oxidizer
                const scalar To() const;


    private:

        //- species
        stringField species_;


        //- mass fraction

            //- mass fraction [-]
            scalarField Y_;

            //- mass fraction in pure fuel [-]
            scalarField Yf_;

            //- mass fraction in pure oxidizer [-]
            scalarField Yo_;


        //- mol fraction

            //- mol fraction [-]
            scalarField X_;

            //- mol fraction in pure fuel [-]
            scalarField Xf_;

            //- mol fraction in pure oxidizer [-]
            scalarField Xo_;


        //- temperature [K]
        scalar T_;

        //- temperature of pure fuel
        scalar Tf_;

        //- temperature of pure oxidizer
        scalar To_;

        //- concentration [-]
        scalar Con_;

        //- is part of fuel mixture
        std::vector<bool> fuel_;

        //- is part of oxidizer mixture
        std::vector<bool> oxidizer_;

        //- Chemistry class reference
        const Chemistry* pChemistry;
};

#endif // MixtureFraction_HPP
