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
    TKC::Properties

Description
    Abstract TKC::Properties class for properties data and calculation

SourceFiles
    properties.cpp

\*---------------------------------------------------------------------------*/

#ifndef Properties_hpp
#define Properties_hpp

#include "propertiesCalc.hpp"
#include "thermo.hpp"
#include "chemistry.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                            Class Properties Declaration
\*---------------------------------------------------------------------------*/

class Properties
:
    public PropertiesCalc
{
    private:

        // Debug
        const bool debug{false};


    public:

        //- Constructor
        Properties(const string, const Thermo&, const Chemistry&);

        //- Destructor
        ~Properties();


        // Time Related Functions

            //- Update the current time stamp [s]
            /*void updateCurrentTime(const scalar);

            //- Return runTime [s]
            scalar runTime() const;

            //- Return deltat
            scalar deltat() const;

            //- Return the current time [s]
            scalar currentTime() const;

            //- Return write control for time
            scalar write() const;
            */


        // Return Functions



        // Summary function

           //- Build the output file that contains all data
           //  void summary(ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Properties_hpp included

// ************************************************************************* //
