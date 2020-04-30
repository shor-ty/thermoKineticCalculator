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
    TKC::Interpreter

Description
    Abstract TKC::Interpreter class that interpretes all data. The class
    inherits just one function which is building the output string and calles
    the summary functions of the classes

SourceFiles
    interpreter.cpp

\*---------------------------------------------------------------------------*/

#ifndef Interpreter_hpp
#define Interpreter_hpp

#include "transport.hpp"
#include "thermo.hpp"
#include "chemistry.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                            Class Numerics Declaration
\*---------------------------------------------------------------------------*/

class Interpreter
{

    public:

        //- Constructor
        Interpreter();

        //- Destructor
        ~Interpreter();


        // Member functions

            //- Create a summary of the loaded data
            void summary(Transport&, const Thermo&, const Chemistry&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Interpreter_hpp included

// ************************************************************************* //
