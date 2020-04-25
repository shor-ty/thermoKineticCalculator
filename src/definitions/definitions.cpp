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

#include "definitions.hpp"

// * * * * * * * * * * * * * * * * Definitions * * * * * * * * * * * * * * * //

std::ostream& TKC::Info = std::cout;

std::ostream& TKC::Error = std::cerr;

std::basic_ostream<char>& (& TKC::endl)(std::basic_ostream<char>&) = std::endl;


void TKC::ErrorMsg
(
    const string msg,
    const char* file,
    const unsigned long line
)
{
    Error<< "\n    *** Error in " << file << " line " << line << "\n"
         << "    " << msg << "\n\n"
         << "    If there is a bug or a problem that you can not solve,\n"
         << "    do not hesitate to write an email to "
         << "Tobias.Holzmann@Holzmann-cfd.de.\n" << endl;

    std::terminate();
}


void TKC::Warning(const string msg, const char* file, const unsigned long line)
{
    Error<< "\n    * Warning in " << file << " line " << line << "\n" << msg
         << "\n" << endl;
}


void TKC::NotImplemented(const char* file, const size_t line)
{
    Error<< "\n"
         << "    * The functionality is not implemented.\n" << endl;
}


TKC::string TKC::Header()
{
    string header = \
"\
/*------------------------------------------------------------------------*\\\
\n|  c-o-o-c-o-o-o             |                                             |\
\n|  |     |     T hermo       | TKC: The Open Source Thermo-Kinetic Library |\
\n|  c-o-o-c     K inetic      | Version: 1.0.0                              |\
\n|  |     |     C alculator   | Web: www.Holzmann-cfd.com                   |\
\n|  c     c-o-o-o             |                                             |\
\n\
\\*------------------------------------------------------------------------*/\
\n";

    return header;
}

// * * * * * * * * * * * * * * String Conversation * * * * * * * * * * * * * //



void TKC::Footer(const scalar startTime)
{
    const scalar execTime = (clock()-startTime) / (scalar) CLOCKS_PER_SEC;

    Info<< "\n\n c-o Execution: " << execTime  << " s\n\n"
        << " ============================================================\n\n"
        << " c-o Programmed by Tobias Holzmann\n\n"
        << " c-o Tobias.Holzmann@Holzmann-cfd.de\n" << endl;
}


// ************************************************************************* //
