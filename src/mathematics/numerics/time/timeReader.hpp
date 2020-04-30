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
    TKC::TimeReader

Description
    Abstract TKC::TimeReader class for building and calculating matrices

SourceFiles
    numerics.cpp

\*---------------------------------------------------------------------------*/

#ifndef TimeReader_hpp
#define TimeReader_hpp

#include "definitions.hpp"
#include "time.hpp"
#include "stringManipulator.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                            Class TimeReader Declaration
\*---------------------------------------------------------------------------*/

class TimeReader
:
    public StringManipulator
{
    private:

        // Private Data

            //- File in which the data are stored
            string file_;


    public:

        // Constructor and Destructor

            //- Constructor
            TimeReader(const string&);

            //- Destructor
            ~TimeReader();


        // Member Functions

            //- Read corresponding file and insert data to the Time class
            void read(Time&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // TimeReader_hpp included

// ************************************************************************* //
