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
    AFC::ChemistryCalc
    
Description
    Abstract AFC::ChemistryCalc class for chemistry calculation

SourceFiles
    chemistry.cpp

\*---------------------------------------------------------------------------*/

#ifndef ChemistryCalc_hpp
#define ChemistryCalc_hpp

#include "chemistryData.hpp"
#include "thermo.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class ChemistryCalc Declaration
\*---------------------------------------------------------------------------*/

class ChemistryCalc
{
    private:

        // Debug 
        bool debug{false};

    public:

        //- Constructor
        ChemistryCalc();

        //- Destructor
        ~ChemistryCalc();


        // Calculation Functions

            //- Calculate the reaction rates kf and kb
            void k
            (
                const scalar&,
                const map<word, scalar>&,
                const Thermo&,
                ChemistryData&
            );

            //- Calculate kf
            void calculateKf
            (
                scalarField&,
                const int&,
                const scalar&,
                const map<word, scalar>&,
                ChemistryData&
            );

            //- Calculate kb
            void calculateKb
            (
                scalarField&,
                const scalarField&,
                const int&,
                const scalar&,
                const map<word, scalar>&,
                const Thermo&,
                ChemistryData&
            );

            //- Calculate [M] partner
            scalar calcM
            (
                const map<word, scalar>&,
                const int&
            );

            //- Calculate [M] partner out of [10]
            scalar calcM_Warnatz
            (
                const map<word, scalar>&
            );

            //- Calculate k with standard arrhenius
            scalar arrhenius
            (
                const scalar&,
                const scalar&,
                const scalar&,
                const scalar&
            );

            //- Calculate the source term of each species (omega)
            void omega
            (
                const map<word, scalar>&,
                ChemistryData&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ChemistryCalc_hpp included

// ************************************************************************* //
