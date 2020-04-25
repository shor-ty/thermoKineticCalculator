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
    TKC::Rosenbrock
    
Description
    Abstract TKC::Rosenbrock class. Numerical implementation based on 
    W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery
    Numerical Recipes - The Art of Scientific Computing - Third Edition

SourceFiles
    numerics.cpp

\*---------------------------------------------------------------------------*/

#ifndef Rosenbrock_hpp
#define Rosenbrock_hpp

#include "definitions.hpp"
#include "matrix.hpp"
#include "vector.hpp"
#include "ODE.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                            Class Rosenbrock Declaration
\*---------------------------------------------------------------------------*/

class Rosenbrock
:
    public ODE
{
    private:

        // Debug switch
        bool debug_{true};

        //- Solutions Vectors 
        //Vector dym_;
        //Vector dyt_;
        //Vector yt_;

        //- Helper quantities

    public:

        //- Constructor 
        Rosenbrock
        (
            Chemistry& 
        );

        //- Destructor
        ~Rosenbrock();


        // Member functions

            //- Solve the ODE
            void solve
            (
                const scalar,
                const wordList&,
                const map<word, scalar>&
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Rosenbrock_hpp included

// ************************************************************************* //
