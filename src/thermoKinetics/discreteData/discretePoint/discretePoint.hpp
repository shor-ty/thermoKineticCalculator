/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     T hermo       | Open Source Thermo-Kinetic Library
  c-o-o-c     K iknetic     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Constructor.

    TKC is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    TKC is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TKC; if not, see <http://www.gnu.org/licenses/>

Class
    TKC::DiscretePoint

Description
    Abstract TKC::DiscretePoint class for a single discrete point. The class
    stores all data relevant thermo-kinetic data for one discrete point.

SourceFiles
    discretePoint.cpp

\*---------------------------------------------------------------------------*/

#ifndef DiscretePoint_hpp
#define DiscretePoint_hpp

#include "discretePoint.hpp"
#include "definitions.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                            Class DiscretePoint Declaration
\*---------------------------------------------------------------------------*/

class DiscretePoint
{
    private:

        // Private Data

            //- Concentration of each species [g/mol]
            map<word, scalar> C_;

            //- Mass fraction of each species [-]
            map<word, scalar> Y_;

            //- Mole fraction of each species [-]
            map<word, scalar> X_;

            //- Temperature at point [K]
            scalar T_;

            //- Mean density [kg/m^3]
            scalar rho_;

            //- Mean molecular weight [g/mol]
            scalar MMW_;



    public:

        //- Constructor
        DiscretePoint();

        //- Destructor
        ~DiscretePoint();


        // Insert Functions

            //- Concentration of species s [g/mol]
            void C(const word, const scalar);

            //- Mass fraction of species s [-]
            void Y(const word, const scalar);

            //- Mole fraction of species s [-]
            void X(const word, const scalar);

            //- Temperature [K]
            void T(const scalar);

            //- Mean density [kg/m^3]
            void rho(const scalar);

            //- Mean molecular weight [g/mol]
            void MMW(const scalar);


        // Return Functions

            //- Concentration of all species [g/mol]
            const map<word, scalar>& C() const;

            //- Concentration of species s [g/mol]
            const scalar& C(const word) const;

            //- Mass fraction of all species [-]
            const map<word, scalar>& Y() const;

            //- Mass fraction of species s [-]
            const scalar& Y(const word) const;

            //- Mole fraction of all species [-]
            const map<word, scalar>& X() const;

            //- Mole fraction of species s [-]
            const scalar& X(const word) const;

            //- Temperature [K]
            const scalar& T() const;

            //- Mean density [kg/m^3]
            const scalar& rho() const;

            //- Mean molecular weight [g/mol]
            const scalar& MMW() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // DiscretePoint_hpp included

// ************************************************************************* //
