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

#include <time.h>
#include <vector>
#include <string>
#include <sstream>
#include <string.h>
#include <iostream>
#include <iomanip> 
#include <memory>
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

class MixtureFraction; 

/*---------------------------------------------------------------------------*\
                          AFC Defintions and typedef
\*---------------------------------------------------------------------------*/

using ostream = std::ostream;

using fstream = std::fstream;

using scalar = long double;

using boolList = std::vector<bool>;

using scalarField = std::vector<scalar>;

using scalarList = std::vector<scalar>;

using matrix = std::vector<scalarField>;

using word = std::string;

using wordList = std::vector<word>;

using intList = std::vector<int>;

using wordMatrix = std::vector<wordList>;

using wordScalarList = std::vector<std::map<word, scalar> >;

using string = std::string;

using stringList = std::vector<string>;

using lookUpTable = std::vector<std::vector<MixtureFraction> >;

extern std::ostream& Info;

extern std::ostream& Error;

extern std::basic_ostream<char>& (&endl)(std::basic_ostream<char>&);

template<class T> using List = std::vector<T>;

template<class T> using vector = std::vector<T>;

template<class T> using smartPtr = std::unique_ptr<T>;

template<class T, class Z> using map = typename std::map<T, Z>;

template<class T, class Z> using mapList = typename std::vector<map<T, Z> >;

#define forEach(Field, i) for (unsigned int i=0; i<Field.size(); i++)
//#define forAll(Field, entry) for (const auto& entry : Field)
#define forAll(Field, entry) for (auto& entry : Field)
//#define forMap(Field, entry) for (const auto& entry : Field)
//#define forMap(Field, entry) for (auto entry : Field)

void Print();

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

string Header();

void Footer
(
    const scalar& 
);

// * * * * * * * * * * * * * * String Conversation * * * * * * * * * * * * * //

template <typename T>
string toStr
(
    const T& tmp
)
{
    return std::to_string(tmp);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // TYPEDEF_HPP_INCLUDED

// ************************************************************************* //
