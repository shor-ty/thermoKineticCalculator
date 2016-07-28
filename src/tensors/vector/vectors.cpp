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
    const size_t rows = 1,
    const size_t cols = 1,
    const scalar value
)
:
    Tensor(rows, cols, value)
{
    if (debug_)
    {
        Info<< "Constructor Vector (row, value)\n" << endl;
    }

    //- Just avoid to make a matrix
    if (rows == cols)
    {
        FatalError
        (
            "    The rows and cols of the vector are identical. It is\n"
            "    better to use a scalar or a matrix object. This depend\n"
            "    on what you need.",
            __FILE__,
            __LINE__
        );
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



// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

AFC::Vector AFC::Vector::T() const
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
}


// * * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * //

//AFC::Vector AFC::Vector

// ************************************************************************* //
