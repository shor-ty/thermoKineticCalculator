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
    AFC::Vector
    
Description
    Abstract AFC::Vector class that inherits all vector calculations 


SourceFiles
    vector.cpp

\*---------------------------------------------------------------------------*/

#ifndef Vector_hpp
#define Vector_hpp

#include "typedef.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Vector Declaration
\*---------------------------------------------------------------------------*/

class Vector
{
    private:

        // Debug swithc
        bool debug_{true};

        // Private Data


    public:

        //- Constructor 
        Vector();

        //- Constructor that creates a row vector n x 1 and init with the 
        //  given value. Cols vectors are also possible with initialisation
        //  like Vector b(1,4) this will make a column vector 1 x 4
        Vector
        (
            const size_t,
            const size_t cols = 1,
            const scalar value = 0
        );

        //- Destructor
        ~Vector();

        // Operator Functions


        // Member Functions

        // Calculation Functions

            //- Transpose the vector b -> b^T
            Vector T() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Vector_hpp included

// ************************************************************************* //
