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
    AFC::Numerics
    
Description
    Abstract AFC::Numerics class for building and calculating matrices

SourceFiles
    numerics.cpp

\*---------------------------------------------------------------------------*/

#ifndef Numerics_hpp
#define Numerics_hpp

#include "definitions.hpp"
#include "mixtureFraction.hpp"
#include "rosenbrock.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Numerics Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Numerics
{
    private:

        // Debug switch
        bool debug_{true};

        //- Pointer to the ODE solver
        Type* ode_;

        //- Initial chemistry time step for each discrete point Z 
        //  This time step will change during the calculation
        scalarField deltaTChem_;


    public:

        //- Constructor 
        Numerics
        (
            Chemistry&  
        );

        //- Constructor usind make_unique with Type default
        // Numerics : Type (std::make_unique<Type>()) {};




        //- Destructor
        ~Numerics();


        // Member functions

            //- Initialize species and temperature; here we just solve the
            //  laplace equation to get the linear profile (this could also
            //  be done in the initialization (would be much faster)
            //  The scalar dissipatrion rate is set automatically to get
            //  a very fast solution
            void solveForInitialSolution
            (
                MixtureFraction&
            );

            //- Solve the transient laminar flamelet equations for the species
            //  and the temperature
            void solveAdiabaticFlamelet 
            (
                MixtureFraction&
            );
            
            //- Solve the chemistry with Seulex algorithm
            scalar solveChemistry
            (
                const scalar,
                MixtureFraction&
            );

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

            //- Solve the problem Ax = b
            void solveMatrix
            (
                matrix&,
                scalarField&,
                scalarField&
            );

            //- Calculate and return the residual
            scalar residual
            (
                const scalarField&,
                const scalarField&
            ) const;


        // Compare functions

            //- Return maximum of the field
            scalar max
            (
                const scalarField&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Numerics_hpp included

// ************************************************************************* //
