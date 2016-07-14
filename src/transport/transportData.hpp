/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     |
  |     |     C onstructor  | Copyright (C) 2015 Holzmann-cfd
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Constructor.

    AFC is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or 
    (at your option) any later version.

    AFC is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AFC; if not, see <http://www.gnu.org/licenses/>

Class
   AFC::TransportData
   
Description
    This class contains all transport data 

SourceFiles
    transportData.cpp

\*---------------------------------------------------------------------------*/

#ifndef TransportData_hpp
#define TransportData_hpp

#include "typedef.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*--------------------------------------------------------------------------*\
                         Class TransportData Declaration
\*--------------------------------------------------------------------------*/

class TransportData
{
    private:

        // Private data

            //- List of species
            wordList species_;

            //- List of species used in the elementar reactions
            wordList chemSpecies_;

            //- Hashtable for geometrical configuration
            //  Definiton:
            //      + 0 molecule is single atom
            //      + 1 molecule is linear
            //      + 2 molecule is nonlinear
            map<word, int> geoConfig_;

            //- Hashtable for Lennard-Jones potential well depth eps/kb [K]
            map<word, scalar> lenJonPot_;

            //- Hashtable for Lennard-Jones collision dimater C in Angstroms
            map<word, scalar> lenJonCollDia_;

            //- Hashtable for dipole moment mu in debey
            //  a debey is 10^-18 [cm^3/2]
            map<word, scalar> dipMom_;

            //- Hashtable for polarizability alpha in cubic Angstroms
            map<word, scalar> pol_;

            //- Hashtable for rotational relaxation collision number Zrot
            //  at 298K
            map<word, scalar> rotRelCollNumb_;


        // Transport data for fitting procedure 

            //- Pure species viscosity
            //  Species, Temp, Value
            //mapMap<word, scalar, scalar> viscosity_;

            //- Pure binary diffusion coefficients


    public:

        //- Constructor 
        TransportData();

        //- Destructor
        ~TransportData();


        // Member functions



        // Insert functions, from TransportReader:: delegated 

            //- Insert species
            void insertSpecies
            (
                const word&
            );

            //- Insert geometrical configuration
            void insertGeoConfig
            (
                const int&
            );

            //- Insert Lennard-Jones potential
            void insertLenJonPot
            (
                const scalar&
            );

            //- Insert Lennard-Jones collision diameter
            void insertLenJonCollDia
            (
                const scalar&
            );

            //- Insert dipole momentum
            void insertDipMom
            (
                const scalar&
            );

            //- Insert polarizability
            void insertPol
            (
                const scalar&
            );

            //- Insert rotational relaxation collision number
            void insertRotRelCollNumb
            (
                const scalar&
            );

        // Insert function from afc.cpp

            //- Insert chemical species
            void chemSpecies
            (
                const wordList&
            );


        // Return functions
        
            //- Return species 
            wordList species() const;

            //- Return chemical species
            wordList chemSpecies() const;

            //- Return the geometrical configuration
            int geometricalConfig
            (
                const word&
            );

            //- Return Lennard-Jones collision diameter [Angstroms]
            scalar LJCD
            (
                const word&
            ) const;

            //- Return Lennard-Jones potential well depth eps/kb [K]
            scalar LJP
            (
                const word&
            ) const;

            //- Return the dipole moment [debey]
            scalar muk
            (
                const word&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // TransportData_hpp included

// ************************************************************************* //
