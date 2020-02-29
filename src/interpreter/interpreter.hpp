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
    AFC::Interpreter
    
Description
    Abstract AFC::Interpreter class that interpretes all data. The class
    inherits just one function which is building the output string and calles
    the summary functions of the classes

SourceFiles
    interpreter.cpp

\*---------------------------------------------------------------------------*/

#ifndef Interpreter_hpp 
#define Interpreter_hpp 

#include "chemistry.hpp"
#include "thermo.hpp"
#include "transport.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Numerics Declaration
\*---------------------------------------------------------------------------*/

class Interpreter 
{
    private:

        //- Switch if interpreter is used or not
        bool interprete_{false};


    public:

        //- Constructor 
        Interpreter();

        //- Destructor
        ~Interpreter();


        // Member functions

            bool analyze() const;

            void analyze(const bool&);

            //- Create a summary of the loaded data
            void summary 
            (
                const Chemistry&,
                const Thermo&,
                const Transport&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Interpreter_hpp included

// ************************************************************************* //
