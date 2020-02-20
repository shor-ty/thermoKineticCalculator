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
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Interpreter::~Interpreter()
{}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //


bool AFC::Interpreter::analyze() const
{
    return interprete_;
}


void AFC::Interpreter::analyze
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
    const Transport& transport
)
{
    //- Create folder  
    system("mkdir -p summary");

    std::filebuf file;

    //- Chemistry summary
    {
        Info<< " c-o Interprete chemistry data\n" << endl;

        file.open("summary/chemistry.afc", std::ios::out);

        std::ostream data(&file);

        //- Set scientific notation
        data.setf(std::ios::scientific, std::ios::floatfield);

        //- Build the chemistry summary
        chemistry.summary(data);

        file.close();
    }

    //- Thermodynamic summary 
    {
        Info<< " c-o Interprete thermodynamic data\n" << endl;

        file.open("summary/thermodynamic.afc", std::ios::out);

        ostream data(&file);

        //- Set scientific notation
        data.setf(std::ios::scientific, std::ios::floatfield);

        //- Build the thermodynamic summary
        thermo.summary(data);

        file.close();
    }

    //- Transport summary 
    {
        Info<< " c-o Interprete transport data\n" << endl;

        file.open("summary/transport.afc", std::ios::out);

        ostream data(&file);

        //- Set scientific notation
        data.setf(std::ios::scientific, std::ios::floatfield);

        //- Build transport summary
        transport.summary(data);

        file.close();
    }

    //- Fiting procedure
    {
        Info<< " c-o Interprete polynomial fit\n" << endl;

        file.open("summary/fittingProcedure.afc", std::ios::out);

        ostream data(&file);

        //- Set scientific notation
        data.setf(std::ios::scientific, std::ios::floatfield);
        
        //- Build transport summary
        transport.summaryFittingProcedure(data);

        file.close();
    }

    Info<< " c-o Interperted all data\n" << endl;
}


// ************************************************************************* //
