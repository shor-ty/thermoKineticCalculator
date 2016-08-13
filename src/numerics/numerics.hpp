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
    AFC::Numerics
    
Description
    Abstract AFC::Numerics class for building and calculating matrices

SourceFiles
    numerics.cpp

\*---------------------------------------------------------------------------*/

#ifndef Numerics_hpp
#define Numerics_hpp

#include "typedef.hpp"
#include "mixtureFraction.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Numerics Declaration
\*---------------------------------------------------------------------------*/

class Numerics
{
    private:

        // Debug switch
        bool debug_{false};


    public:

        //- Constructor 
        Numerics();

        //- Destructor
        ~Numerics();


        // Member functions

            //- Solve laplace term by the Finite Difference Method
            //  using a 2nd order discretization scheme for for non-equal
            //  spaced 1D problems [Holzmann] 
            scalar FDMLapacian2ndOrder
            (
                const scalar&,
                const scalar&,
                const scalar&,
                const scalar&,
                const scalar&,
                const scalar&
            ) const;

            //- Solve the Flamelet equations
            void solveFlamelet
            (
                MixtureFraction&,
                const scalar&,
                const scalar&
            );

            //- solve Jacobian
            void jacobian
            (
                MixtureFraction&
            );

            //- Solve the problem Ax = b
            void solveMatrix
            (
                matrix&,
                scalarField&,
                scalarField&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Numerics_hpp included

// ************************************************************************* //
