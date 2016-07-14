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

\*---------------------------------------------------------------------------*/

#include <fstream>
#include <iomanip>
#include "interpreter.hpp"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

AFC::Interpreter::Interpreter()
{
    if (debug)
    {
        Info<< "Constructor Interpreter\n" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Interpreter::~Interpreter()
{
    if (debug)
    {
        Info<< "Destructor Interpreter \n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //


bool AFC::Interpreter::analyse() const
{
    return interprete_;
}


void AFC::Interpreter::analyse
(
    const bool& interprete 
)
{
    interprete_ = interprete;
}


void AFC::Interpreter::summary
(
    const Chemistry& chemistry,   
    const Thermo& thermo,
    const Transport& Transport
)
{
    //- Create folder  
    system("mkdir -p summary");

    std::filebuf file;

    //- Chemistry summary
    {
        file.open("summary/chemistry.afc", std::ios::out);

        std::ostream data(&file);

        //- Set scientific notation
        data.setf(std::ios::scientific, std::ios::floatfield);

        //- Header
        data<< Header() << "\n"; 

        //- Build the chemistry summary
        chemistry.summary(data);

        file.close();
    }

    //- Thermodynamic summary 
    {
        file.open("summary/thermodynamic.afc", std::ios::out);

        std::ostream data(&file);

        //- Set scientific notation
        data.setf(std::ios::scientific, std::ios::floatfield);

        //- Header
        data<< Header() << "\n"; 

        //- Build the thermodynamic summary
        thermo.summary(data);

        file.close();
    }

    //- Transport summary 
    {
        file.open("summary/transport.afc", std::ios::out);

        std::ostream data(&file);

        //- Set scientific notation
        data.setf(std::ios::scientific, std::ios::floatfield);

        //- Header
        data<< Header() << "\n"; 

        file.close();
    }

    Info<< " c-o Interperted all data\n" << endl;
}


// ************************************************************************* //
