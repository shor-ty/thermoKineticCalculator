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
    AFC::Transport
    
Description
    Abstract AFC::Transport class for transport data and calculation

SourceFiles
    transport.cpp

\*---------------------------------------------------------------------------*/

#ifndef Transport_hpp
#define Transport_hpp

#include "transportCalc.hpp"
#include "thermo.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Transport Declaration
\*---------------------------------------------------------------------------*/

class Transport
:
    public TransportCalc
{
    private:

        // Private reference data


    public:

        //- Constructor
        Transport(const string, const Thermo&);

        //- Destructor
        ~Transport();


        // Calculation Functions

            //- Calculate fitting coefficients for polynomials
            void fitCurves();


        // Summary functions

            //- Build the summary
            void summary(ostream&) const;

            //- Build the fitting procedure summary
            void summaryFittingProcedure(ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Transport_hpp included

// ************************************************************************* //
