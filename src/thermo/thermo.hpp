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
    AFC::Thermo
    
Description
    Abstract AFC::Thermo class for thermo data and calculation

SourceFiles
    thermo.cpp

\*---------------------------------------------------------------------------*/

#ifndef Thermo_hpp
#define Thermo_hpp

#include "thermoReader.hpp"
#include "thermoData.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Thermo Declaration
\*---------------------------------------------------------------------------*/

class Thermo
{
    private:

        // Private pointer data

            //- ThermoData object
            ThermoData thermoData_;
    

    public:

        //- Constructor with fileName
        Thermo
        (
            const string&,
            const bool&
        );

        //- Destructor
        ~Thermo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Thermo_hpp included

// ************************************************************************* //
