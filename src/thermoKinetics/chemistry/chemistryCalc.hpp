/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

class Thermo;

/*---------------------------------------------------------------------------*\
                            Class ChemistryCalc Declaration
\*---------------------------------------------------------------------------*/

class ChemistryCalc
:
    public ChemistryData
{

    private:

        // Private reference data

            //- Reference tp the thermo object
            const Thermo& thermo_;


    public:

        //- Constructor
        ChemistryCalc(const string, const Thermo&);

        //- Destructor
        ~ChemistryCalc();


        // Member Functions

            //- Return the reference to the Thermo object
            const Thermo thermo() const;


        // Calculation Functions

            //- Calculate the reaction rates kf and kb [units depend]
            void kfkb
            (
                const int,
                const scalar,
                const map<word, scalar>&
            );

            //- Calculate reaction rate kf
            scalar kf
            (
                const int,
                const scalar,
                const map<word, scalar>& = map<word, scalar>(),
                const bool kb = false
            ) const;

            //- Calculate reaction rate kb
            scalar kb
            (
                const int,
                const scalar,
                const map<word, scalar>& = map<word, scalar>()
            ) const;

            //- Calculate equilibrium reaction rate keq
            scalar keq(const int, const scalar) const;

            //- Calculate k with standard arrhenius [units depend]
            scalar arrhenius
            (
                const scalar,
                const scalar,
                const scalar,
                const scalar
            ) const;

            //- Lindemann approach
            scalar Lindemann
            (
                const int,
                const scalar,
                const map<word, scalar>&
            ) const;

            //- Calculate Fcent for TROE formulation
            scalar Fcent(const int, const scalar) const;

            //- Calculate Flog for TROE formulation
            scalar Flog(const int, const scalar, const scalar) const;

            //- Calculate the source term of species s (omega)
            scalar omega
            (
	            const word,
                const scalar,
                const map<word, scalar>&
            ) const;

            //- Calculate [M] partner [mol/m^3]
            scalar M(const int, const map<word, scalar>&) const;

            //- Calculate dH for reaction r and given temperature
            scalar dh(const int, const scalar) const;

            //- Calculate dG for reaction r and given temperature
            scalar dg(const int, const scalar) const;

            //- Calculate dS for reaction r and given temperature
            scalar ds(const int, const scalar) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ChemistryCalc_hpp included

// ************************************************************************* //
