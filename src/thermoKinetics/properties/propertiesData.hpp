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
    TKC::PropertiesData

Description
    Abstract TKC::PropertiesData class for properties data and calculation

SourceFiles
    properties.cpp

\*---------------------------------------------------------------------------*/

#ifndef PropertiesData_hpp
#define PropertiesData_hpp

#include "definitions.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                            Class PropertiesData Declaration
\*---------------------------------------------------------------------------*/

class PropertiesData
{
    private:

        // Debug
        const bool debug{false};


        // General private data

            //- Inert species
            word inert_{"N2"};

            //- Initial composition n [g/mol]
            map<word, scalar> n_;

            //- Initial mol fraction X [-]
            map<word, scalar> X_;

            //- Initial mass fraction Y [-]
            map<word, scalar> Y_;

            //- initial temperature [K]
            scalar T_{0};

            //- Pressure at which the calculation takes place [Pa]
            scalar p_{0};

            
        // Boolean

            //- Input either mole or mass fraction or concentration
            bool inputMole_{false};
            bool inputMass_{false};
            bool inputConcentration_{false};


    public:

        //- Constructor
        PropertiesData(const string);

        //- Destructor
        ~PropertiesData();


        // Member functions

            //- Set inert species
            void insertInertSpecies(const word);

            //- Insert initial concentration composition [g/mol]
            void insertC(const word, const scalar);

            //- Insert initial mole composition [-]
            void insertX(const word, const scalar);

            //- Insert initial mass fraction [-]
            void insertY(const word, const scalar);

            //- Insert initial temperature [K] 
            void insertT(const scalar);

            //- Insert pressure [Pa]
            void insertP(const scalar);

            //- Insert input mode (mole, mass, concentration)
            void inputMode(const word);


        // Return Functions

            //- Return inert species
            const word inertSpecies() const;

            //- Return initial composition concentration [g/mol]
            const map<word, scalar> C() const;

            //- Return initial mole fraction
            const map<word, scalar> X() const;

            //- Return initial mass fraction
            const map<word, scalar> Y() const;

            //- Return initial temperature [K]
            const scalar T() const;

            //- Return pressure [Pa]
            const scalar p() const;

            //- Return word of input mode (mole, mass or concentration)
            const word inputMode() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // PropertiesData_hpp included

// ************************************************************************* //
