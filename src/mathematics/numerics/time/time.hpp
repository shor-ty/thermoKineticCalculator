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
    TKC::Time

Description
    Abstract TKC::Time class for building and calculating matrices

SourceFiles
    numerics.cpp

\*---------------------------------------------------------------------------*/

#ifndef Time_hpp
#define Time_hpp

#include "definitions.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                            Class Time Declaration
\*---------------------------------------------------------------------------*/

class Time
{
    private:

        //- Time related data

            //- Actual run time
            scalar runTime_{0};

            //- End time
            scalar endTime_{0};

            //- Write results
            scalar writeTime_{0};

            //- Kinetic time step
            scalar deltaTKinetic_{0};

            //- Flow time step (for Flamelet calculation)
            scalar deltaTFlow_{0};


    public:

        // Constructors and Destructors

            //- Constructor
            Time(const string);

            //- Destructor
            ~Time();


        // Operators

            //- Increment runTime based on a given delta t
            void operator+=(const scalar);

        // Insert Functions

            //- Insert the end time for calculation
            void endTime(const scalar);

            //- Insert the write time for the results
            void writeTime(const scalar);

            //- Update the kinetic time step
            void dTKinetic(const scalar);

            //- Update the flow time step
            void dTFlow(const scalar);


        // Return Functions

            //- Return the actual run time
            const scalar runTime() const;

            //- Return the end time of the calculation
            const scalar endTime() const;

            //- Return the time after which the results are stored
            const scalar writeTime() const;

            //- Return the actual time step for the kinetic
            const scalar dTKinetic() const;

            //- Return the actual time step for the flow
            const scalar dTFlow() const;


        // Member functions

            //- Returns true if the end time is not reached
            const bool loop() const;

            //- Check the data
            void checkData() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Time_hpp included

// ************************************************************************* //
