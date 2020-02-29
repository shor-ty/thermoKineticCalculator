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

\*---------------------------------------------------------------------------*/

#include "typedef.hpp"

// * * * * * * * * * * * * * * * * Definitions * * * * * * * * * * * * * * * //

std::ostream& AFC::Info = std::cout;

std::ostream& AFC::Error = std::cerr;

std::basic_ostream<char>& (& AFC::endl)(std::basic_ostream<char>&) = std::endl;


void AFC::ErrorMsg
(
    const string msg,
    const char* file, 
    const unsigned long line
)
{
    Error<< "\n    *** Error in " << file << " line " << line << "\n"
         << "    " << msg << "\n\n"
         << "    If there is a bug or a problem you can not solve,\n"
         << "    do not hesitate to write an email to "
         << " Tobias.Holzmann@holzmann-cfd.com.\n" << endl;

    std::terminate();
}


void AFC::Warning(const string msg, const char* file, const unsigned long line)
{
    Error<< "\n    * Warning in " << file << " line " << line << "\n" << msg
         << "\n" << endl;
}


void AFC::NotImplemented(const char* file, const size_t line)
{
    Error<< "\n"
         << "    * The functionality is not implemented.\n" << endl;
}


AFC::string AFC::Header()
{
    string header = \
"\
/*------------------------------------------------------------------------*\\\
\n|  c-o-o-c-o-o-o             |                                             |\
\n|  |     |     A utomatic    | AFC: The Open Source Flamelet Toolbox       |\
\n|  c-o-o-c     F lamelet     | Version: 1.0.0                              |\
\n|  |     |     C onstructor  | Web: www.Holzmann-cfd.com                   |\
\n|  c     c-o-o-o             |                                             |\
\n\
\\*------------------------------------------------------------------------*/\
\n";

    return header;
}

// * * * * * * * * * * * * * * String Conversation * * * * * * * * * * * * * //



void AFC::Footer(const scalar startTime)
{
    const scalar execTime = (clock()-startTime) / (scalar) CLOCKS_PER_SEC; 

    Info<< "\n\n c-o Execution: " << execTime  << " s\n\n" 
        << " ============================================================\n\n"
        << " c-o Programmed by Tobias Holzmann\n\n"
        << " c-o Tobias.Holzmann@Holzmann-cfd.de\n" << endl;
}


// ************************************************************************* //
