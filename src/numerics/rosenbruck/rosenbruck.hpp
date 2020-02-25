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
    AFC::Rosenbruck
    
Description
    Abstract AFC::Rosenbruck class. Numerical implementation based on 
    W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery
    Numerical Recipes - The Art of Scientific Computing - Third Edition

SourceFiles
    numerics.cpp

\*---------------------------------------------------------------------------*/

#ifndef Rosenbruck_hpp
#define Rosenbruck_hpp

#include "typedef.hpp"
#include "matrix.hpp"
#include "jacobian.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Rosenbruck Declaration
\*---------------------------------------------------------------------------*/

class Rosenbruck
{
    private:

        // Debug switch
        bool debug_{false};


    public:

        //- Constructor 
        Rosenbruck();

        //- Destructor
        ~Rosenbruck();


        // Member functions

            //- Solve using euler algorithm
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Rosenbruck_hpp included

// ************************************************************************* //
