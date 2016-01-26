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

            void k
            (
                const ChemistryData&
            ) const;

            //- Calculate the reaction rates kf and kb [units depend]
            void calculatekfkb 
            (
                const int&,
                scalar&,
                scalar&,
                const scalar&,
                const map<word, scalar>&,
                const Thermo&,
                const ChemistryData&
            );

            //- Calculate kf
            scalar calculateKf 
            (
                const int&,
                const scalar&,
                const scalar&,
                const ChemistryData&
            );

            //- Calculate kb
            scalar calculateKb 
            (
                const scalar&,
                const scalar&
            );

            //- Calculate Kc
            scalar calculateKc 
            (
                const int&,
                const scalar&,
                const Thermo&,
                const ChemistryData&
            );

            //- Calculate k with standard arrhenius [units depend]
            scalar arrhenius
            (
                const scalar&,
                const scalar&,
                const scalar&,
                const scalar&
            ) const;

            //- Calculate the source term of each species (omega)
            void calculateOmega
            (
                const scalar&,
                map<word, scalar>&,
                const Thermo&,
                const ChemistryData&
            );
            
            //- Calculate [M] partner [mol/m^3]
            scalar calculateM
            (
                const int&,
                const map<word, scalar>&,
                const ChemistryData&
            );

            //- Get information about TBR (if [M] is included or not) and
            //  if we have to use enhanced factors for calculating [M]
            bool thirdBodyReaction
            (
                const int&,
                bool&,
                const ChemistryData&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ChemistryCalc_hpp included

// ************************************************************************* //
