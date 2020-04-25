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

\*---------------------------------------------------------------------------*/

#include <fstream>
#include <iomanip>
#include "interpreter.hpp"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

TKC::Interpreter::Interpreter()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::Interpreter::~Interpreter()
{}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //


bool TKC::Interpreter::analyze() const
{
    return interprete_;
}


void TKC::Interpreter::analyze
(
    const bool& interprete 
)
{
    interprete_ = interprete;
}


void TKC::Interpreter::summary
(
    const Chemistry& chemistry,   
    const Thermo& thermo,
    const Transport& transport,
    const Properties& properties
)
{
    //- Create folder  
    system("mkdir -p summary");

    std::filebuf file;

    /*
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

    //- Fitting procedure
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

    //- Properties summary 
    {
        Info<< " c-o Interprete chemistry setup (afcDict)\n" << endl;

        file.open("summary/properties.afc", std::ios::out);

        ostream data(&file);

        //- Set scientific notation
        data.setf(std::ios::scientific, std::ios::floatfield);

        //- Build properties summary
        properties.summary(data);

        file.close();
    }

    Info<< " c-o Interperted all data\n" << endl;
    */
}


// ************************************************************************* //
