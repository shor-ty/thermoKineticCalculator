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

Description
    Include definition, makros, types etc.

\*---------------------------------------------------------------------------*/

#ifndef TYPEDEF_HPP
#define TYPEDEF_HPP

#include <vector>
#include <string>
#include <sstream>
#include <string.h>
#include <iostream>
#include <memory>
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

class MixtureFraction; 

/*---------------------------------------------------------------------------*\
                          AFC Defintions and typedef
\*---------------------------------------------------------------------------*/

using scalar = double;

using boolList = std::vector<bool>;

using scalarField = std::vector<double>;

using matrix = std::vector<scalarField>;

using word = std::string;

using wordList = std::vector<word>;

using wordMatrix = std::vector<wordList>;

using string = std::string;

using stringList = std::vector<string>;

using lookUpTable = std::vector<std::vector<std::vector<MixtureFraction> > >;

extern std::ostream& Info;

extern std::ostream& Error;

extern std::basic_ostream<char>& (&endl)(std::basic_ostream<char>&);

template<class T> using vector = std::vector<T>;

template<class T> using smartPtr = std::unique_ptr<T>;

template<class T, class Z> using map = typename std::map<T, Z>;

#define forAll(scalarField, i) for \
    (unsigned int i=0; i<(scalarField).size(); i++)

void FatalError
(
    const string,
    const char*, 
    const unsigned long
);

void Warning
(
    const string,
    const char*, 
    const unsigned long
);

void Header();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // TYPEDEF_HPP_INCLUDED

// ************************************************************************* //
