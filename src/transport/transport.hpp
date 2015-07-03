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

//- user def. headers
#include "../definitions/typedef.hpp"
#include "../database/elements.hpp"


#ifndef TRANSPORT_HPP
#define TRANSPORT_HPP


class Transport
:
    public Elements
{
    public:

        //- constructor
        Transport();

        //- destructor
        virtual ~Transport();


    public:

        //- functions

            //- increment the vectors
            void transportDataIncrement();

            //- read transport file
            void readTransportFile
            (
                const normalString&,
                const stringField&
            );

            //- open file
            stringField openFile
            (
                const normalString&
            );

            //- split string; delimiter is whitespace
            stringField splitString
            (
                const normalString&
            );

            //- split string; delimiter is given
            stringField splitString
            (
                const normalString&,
                const char
            );

            const bool TRANS
            (
                const unsigned int&
            ) const;


    private:

        //- Lennard-Jones potential epsilon/kb in [K]
        scalarField LJP_;

        //- Lennard-Jones collisions diameter sigma [Angstrom]
        scalarField LJCD_;

        //- Dipole moment mu in Debye
        scalarField DPM_;

        //- plorizability alpha in cubic Angstrom
        scalarField P_;

        //- rotational relaxation collision number Zrot at 298K
        scalarField RRCN_;

        std::vector<bool> TRANS_;

};

#endif // TRANSPORT_HPP
