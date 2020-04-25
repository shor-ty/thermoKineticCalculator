/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     T hermo       | Open Source Thermo-Kinetic Library
  c-o-o-c     K iknetic     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Creator.

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

Description
    Example how to use the TKC for solving matrix equations like Ax=b


\*---------------------------------------------------------------------------*/

#include "definitions.hpp" 
#include "matrix.hpp"
#include "vector.hpp"
#include "LUDecompose.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace TKC;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
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

    Info<< "The matrix of the LU decomposition is:\n";
    LUD();

    Info<< "Solving Ax = b using the LU decomposition matrix:\n";
 
    //- Solve Ax = b using LU decomposition
    LUD.solve(b, x);
    x();

    //- Improve the solution
    LUD.improveSolution(b, x);

    Info<< "The solution of the LU decomposition with improved solution:\n";
    x();


    //- Check the solution
    const Vector bb = A * x;

    Info<< "The solution of Ax is equal to (should be vector b):\n";
    bb();


    Info<< "The error of bb - b is:\n" << endl;

    //- bb is a const vector so we 
    const Vector residual = bb - b;
    residual();

    Footer(startTime);

    return 0;
}


// ************************************************************************* //
