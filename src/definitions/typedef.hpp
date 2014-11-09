/*---------------------------------------------------------------------------*\
| =========                |                                                  |
| \\      /  A utomatic    | Holzmann-cfd                                     |
|  \\    /   F lamelet     | Version: 1.0                                     |
|   \\  /    C onstructor  | Web: www.Holzmann-cfd.de                         |
|    \\/                   |                                                  |
\*---------------------------------------------------------------------------*/
/*
»
»
»
»
»
»
»
»
»
»
\*---------------------------------------------------------------------------*/
//- system headers
#include <vector>
#include <string>
#include <iostream>
//- user def. headers

#ifndef TYPEDEF_HPP_INCLUDED
#define TYPEDEF_HPP_INCLUDED

//- scalar field as double
typedef std::vector<double> scalarField;
typedef std::vector<std::string> stringField;
typedef std::string normalString;

//- gas constant [J/molK]
const double R = 8.314461175;

//- Loops
#define forAll(scalarField, i) for (unsigned int i=0; i<(scalarField).size(); i++)

#endif // TYPEDEF_HPP_INCLUDED




