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
    TKC::Thermo

Description
    Abstract TKC::Thermo class for thermo data and calculation

SourceFiles
    thermo.cpp

\*---------------------------------------------------------------------------*/

#ifndef Thermo_hpp
#define Thermo_hpp

#include "thermoCalc.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                            Class Thermo Declaration
\*---------------------------------------------------------------------------*/

class Thermo
:
    public ThermoCalc
{
    private:

        // Private pointer data

            //- Debug switch
            bool debug_{false};


    public:

        //- Constructor with fileName
        Thermo(const string, const bool thermoInChemistry = false);

        //- Destructor
        ~Thermo();


        // Summary function

            //- Build the output file that contains all data
            void summary(ostream&) const;

            //- Build NASA Coefficient table
            void NASAPolynomials(ostream&, const word) const;

            //- Build thermoanalyse table (calcualte thermo properties)
            void thermoTable(ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Thermo_hpp included

// ************************************************************************* //
