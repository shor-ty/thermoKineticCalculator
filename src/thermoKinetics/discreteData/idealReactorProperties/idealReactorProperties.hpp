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
    TKC::IdealReactorProperties

Description
    Abstract TKC::IdealReactorProperties class for properties data and calculation

SourceFiles
    properties.cpp

\*---------------------------------------------------------------------------*/

#ifndef IdealReactorProperties_hpp
#define IdealReactorProperties_hpp

#include "discretePoint.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                            Class IdealReactorProperties Declaration
\*---------------------------------------------------------------------------*/

class IdealReactorProperties
:
    public DiscretePoint
{
    private:

        // Debug
        const bool debug{false};


        // General private data

            //- Inert species
            word inertSpecies_{"N2"};

            //- Pressure at which the calculation takes place [Pa]
            scalar p_{0};

            //- File for thermo data
            word fileThermo_;

            //- File for chemistry data
            word fileChemistry_;

            //- File for transport data
            word fileTransport_;


        // Boolean

            //- Input either mole or mass fraction or concentration
            bool inputMole_{false};
            bool inputMass_{false};
            bool inputConcentration_{false};


    public:

        //- Constructor with file name (path)
        IdealReactorProperties(const string);

        //- Destructor
        ~IdealReactorProperties();


        // Insert functions

            //- Set inert species
            void inertSpecies(const word);

            //- Insert pressure [Pa]
            void p(const scalar);

            //- Insert input mode (mole, mass, concentration)
            void inputMode(const word);

            //- Insert file (path) for thermo file
            void thermo(const word);

            //- Insert file (path) for chemistry file
            void chemistry(const word);

            //- Insert file (path) for transport file
            void transport(const word);


        // Return Functions

            //- Return inert species
            const word inertSpecies() const;

            //- Return pressure [Pa]
            const scalar p() const;

            //- Return word of input mode (mole, mass or concentration)
            const word inputMode() const;

            //- Return file (path) for thermo file
            const word thermo() const;

            //- Return file (path) for chemistry file
            const word chemistry() const;

            //- Return file (path) for transport file
            const word transport() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // IdealReactorProperties_hpp included

// ************************************************************************* //
