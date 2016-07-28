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

\*---------------------------------------------------------------------------*/

#include "typedef.hpp"
#include "vector.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Vector::Vector()
:
    Tensor(0, 0)
{
    if (debug_)
    {
        Info<< "Constructor Vector() \n" << endl;
    }
}


AFC::Vector::Vector
(
    const Vector& vec
)
:
    Tensor(vec.rows(), vec.cols(), vec.values())
{
    if (debug_)
    {
        Info<< "Constructor Vector(const Vector&) Copy Const\n" << endl;
    }
}


AFC::Vector::Vector
(
    const size_t row,
    const scalar value
)
:
    Tensor(row, 1, value)
{
    if (debug_)
    {
        Info<< "Constructor Vector (row, value)\n" << endl;
    }
}


AFC::Vector::Vector
(
    const scalarField& sF
)
:
    Tensor(sF.size(), 1, sF)
{
    if (debug_)
    {
        Info<< "Constructor Vector (row, scalarField)\n" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Vector::~Vector()
{
    if (debug_)
    {
        Info<< "Destructor Vector \n" << endl;
    }
}


// * * * * * * * * * * * * * * Operator Functions  * * * * * * * * * * * * * //

/*AFC::scalar AFC::Vector::operator()
(
    const size_t i
)
{
    //- Row or col vector
    if (rows() > cols())
    {
        return this->AFC::Tensor::operator()(i,0);
    }
    else
    {
        return this->AFC::Tensor::operator()(0,i);
    }
}


void AFC::Vector::operator()
(
    const size_t i,
    const scalar& value
)
{
    //- Row or col vector
    if (rows() > cols())
    {
        *this->operator()(i, 0, value);
    }
    else
    {
        *this->AFC::Tensor::operator()(0, i, value);
    }
}*/


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

/*AFC::Vector AFC::Vector::T() const
{
    //- Rows and cols of the vector
    const size_t& row = this.rows();
    const size_t& col = this.cols();

    //- Make a new vector
    Vector tmp;
    size_t n{0};

    //- If row vector
    if (row > col)
    {
        tmp.resize(1, rows);
        n = row;
    }
    //- col vector
    else
    {
       tmp.resize(cols, 1); 
       n = col;
    }

    //- Transpone the vector 
    for (size_t i = 0; i < n; ++i)
    {
        //- Transpone a row vector to a col vector
        if (row > col)
        {
//            tmp(1, i, this(r) 
        }
        //- Transpone a col vector to a row vector
        else
        {
            
        }
    }

    return tmp;
}*/


// * * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * //

//AFC::Vector AFC::Vector

// ************************************************************************* //
