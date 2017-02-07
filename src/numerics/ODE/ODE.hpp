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
    AFC::ODE
    
Description
    Abstract AFC::ODE class for numeric calculations for the ODE system

SourceFiles
    ODE.cpp

\*---------------------------------------------------------------------------*/

#ifndef ODE_hpp
#define ODE_hpp

#include "typedef.hpp"
#include "chemistry.hpp"
#include "seulex.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class ODE Declaration
\*---------------------------------------------------------------------------*/

class ODE
{
    private:

        //- Debug switch
        bool debug_{false};

        //  TODO make Readable
        //- Maximum allowed iterations that are used for the integration in
        //  the chemistry calculation
        const size_t maxIterations_{10000};

        //- Absolute tolerance 
        const scalar absoluteTolerance_{1e-15};

        //- Relative tolerance 
        const scalar relativeTolerance_{1e-4};
        
        //- Pointer to solver (Seulex)
        Seulex* solver_;


    public:

        //- Constructor 
        ODE
        (
            const Chemistry&
        );

        //- Destructor
        ~ODE();


        // Member functions

            //- Solve the chemistry using Seulex
            void solve
            (
                const scalar,
                const scalar,
                map<word, scalar>&,
                const scalar,
                scalar&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ODE_hpp included

// ************************************************************************* //
