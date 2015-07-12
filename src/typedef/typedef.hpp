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
#include <string.h>
#include <iostream>
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                          AFC Defintions and typedef
\*---------------------------------------------------------------------------*/

typedef double scalar;

typedef std::vector<bool> boolList;

typedef std::vector<double> scalarField;

typedef std::vector<scalarField> matrix;

typedef std::string word;

typedef std::vector<word> wordList;

typedef std::vector<wordList> wordMatrix;

typedef std::string string;

typedef std::vector<string> stringList;

extern std::ostream& Info;

extern std::ostream& Error;

extern std::basic_ostream<char>& (&endl)(std::basic_ostream<char>&);

template<class T> using smartPtr = std::unique_ptr<T>;

#define forAll(scalarField, i) for \
    (unsigned int i=0; i<(scalarField).size(); i++)

void FatalError
(
    const string,
    const char*, 
    const unsigned long
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // TYPEDEF_HPP_INCLUDED

// ************************************************************************* //
