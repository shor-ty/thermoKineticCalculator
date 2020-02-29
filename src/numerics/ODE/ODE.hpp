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
    AFC::ODE
    
Description
    Abstract AFC::ODE class for numeric calculations for the ODE system
    This class is the base and stores relevant numerical data used for 
    the numerical integrator such as Euler, Rosenbrock or SEULEX.

    The ODE class holds a reference to the chemistry object which includes
    all data regarding reactins, concentrations, mass fraction, and so on

    Additionally, the base class includes the jacobian calculation that is
    needed for different equations

SourceFiles
    ODE.cpp

\*---------------------------------------------------------------------------*/

#ifndef ODE_hpp
#define ODE_hpp

#include "typedef.hpp"
#include "stepStatus.hpp"
#include "chemistry.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class ODE Declaration
\*---------------------------------------------------------------------------*/

class ODE
:
    public StepStatus
{
    public:

        //- Reference to the chemistry object
        Chemistry& chem_;

    private:
        /*
        //  TODO make Readable
        //- Maximum allowed iterations that are used for the integration in
        //  the chemistry calculation
        const size_t maxIterations_{10000};

        //- Absolute tolerance 
        const scalar absoluteTolerance_{1e-15};

        //- Relative tolerance 
        const scalar relativeTolerance_{1e-4};
        
        //- Initial chemical time step 
        scalar timeStep_{1e-7};

        //- Time derivative 
        map<word, scalar> dcdt_;

        //- Field for concentration
        map<word, scalar> c_;
        */

        //- Jacobian

         
        //- Adaptive solver
        scalar safeScale_{0.9};
        scalar alphaDec_{0.25};
        scalar minScale_{0.2};


    public:

        //- Constructor 
        ODE
        (
            Chemistry&
        );

        //- Destructor
        ~ODE();


        // Member functions

            //- Calculate the time derivative dc/dt
            void derivative
            (
                const scalar,
                const scalar,
                const map<word, scalar>& 
            );

            //- Calculate the Jacobian matrix
            void jacobian 
            (
                const scalar,
                const scalar,
                const map<word, scalar>& 
            );
           
            //- Derive the elementar reaction based on the species s and return
            //  the value of the derivation
            scalar derivationOfReaction
            (
                const word,
                const word,
                const wordList&,
                const wordList&,
                const map<word, int>&,
                const map<word, int>&,
                const scalar,
                const scalar,
                const map<word, scalar>&
            ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ODE_hpp included

// ************************************************************************* //
