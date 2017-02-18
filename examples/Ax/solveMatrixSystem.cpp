/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2015 Holzmann-cfd
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Creator.

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
    Example how to use the AFC for solving matrix equations like Ax=b


\*---------------------------------------------------------------------------*/

#include "typedef.hpp" 
#include "matrix.hpp"
#include "vector.hpp"
#include "LUDecompose.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace AFC;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main
(
    int argc,
    char** argv
)
{

    const std::clock_t startTime = clock();

    Info<< Header() << endl;


    //- Generate 3x3 matrix
    Matrix A(3,3);

    //- Build matrix
    A(0,0) = 0.000003;
    A(0,1) = 0.213472;
    A(0,2) = 0.332147;

    A(1,0) = 0.215512;
    A(1,1) = 0.375623;
    A(1,2) = 0.476625;

    A(2,0) = 0.173257;
    A(2,1) = 0.663257;
    A(2,2) = 0.625675;

    //- Generate 3x1 matrix (vector)
    Vector b(3);

    //- Build vector
    b(0) = 0.235262;
    b(1) = 0.127653;
    b(2) = 0.285321;

    Info<< "The coefficient matrix A: \n";

    A();

    Info<< "The source vetor b: \n";
    
    b();

    //- Solution vector
    Vector x(3);

    //- Create a new object and decompose A into LU
    LUDecompose LUD(A);

    //- Solve Ax = b using LU decomposition
    LUD.solve(b, x);

    Info<< "The solution of the LU decomposition method:\n";

    x();

    //- Improve the solution
    LUD.improveSolution(b, x);

    Info<< "The solution of the LU decomposition with improved solution:\n";

    x();


    Footer(startTime);

    return 0;
}


// ************************************************************************* //
