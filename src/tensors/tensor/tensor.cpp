/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2015 Holzmann-cfd c     c-o-o-o             |
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
#include "tensor.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Tensor::Tensor()
:
    nRows_(0),
    nCols_(0)
{
    if (debug_)
    {
        Info<< "Constructor Tensor() \n" << endl;
    }

}


AFC::Tensor::Tensor
(
    const size_t rows,
    const size_t cols,
    const scalar value

)
:
    nRows_(rows),
    nCols_(cols),
    values_(rows * cols, value)
{
    if (debug_)
    {
        Info<< "Constructor Tensor (rows, cols, value)\n" << endl;
    }
}


AFC::Tensor::Tensor
(
    const size_t rows,
    const size_t cols,
    const scalarField& sF
)
:
    nRows_(rows),
    nCols_(cols),
    values_(sF)
{
    if (debug_)
    {
        Info<< "Constructor Tensor (rows, cols, scalarField)\n" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Tensor::~Tensor()
{
    if (debug_)
    {
        Info<< "Destructor Tensor\n" << endl;
    }
}


// * * * * * * * * * * * * * * Operator Functions  * * * * * * * * * * * * * //

AFC::scalar AFC::Tensor::operator()
(
    size_t rows,
    size_t cols
) const
{
    return values_[rows * nCols_ + cols];
}


void AFC::Tensor::operator()
(
    size_t rows,
    size_t cols,
    scalar value
)
{
    values_[rows * nCols_ + cols] = value;
}


void AFC::Tensor::operator()() const
{
    Info<< "size: " << nRows_ << endl;
    Info<< "\n";
    for (size_t i = 0; i < nRows_; i++)
    {
        Info<< "| ";
        for (size_t j = 0; j < nCols_; j++)
        {
            Info.precision(2);
            Info<< std::fixed << std::setw(10)
                << this->operator()(i, j);
        }
        Info << " |\n";
    }
    Info<< "\n";
}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

size_t AFC::Tensor::cols() const
{
    return nCols_;
}


size_t AFC::Tensor::rows() const
{
    return nRows_;
}


AFC::scalarField AFC::Tensor::values() const
{
    return values_;
}


// * * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * //


// ************************************************************************* //
