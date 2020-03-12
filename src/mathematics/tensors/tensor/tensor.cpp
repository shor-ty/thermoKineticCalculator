/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD 
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

#include "definitions.hpp"
#include "tensor.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Tensor::Tensor()
:
    nRows_(0),
    nCols_(0)
{}


AFC::Tensor::Tensor(const size_t rows, const size_t cols, const scalar value)
:
    nRows_(rows),
    nCols_(cols),
    values_(rows * cols, value)
{}


AFC::Tensor::Tensor(const size_t rows, const size_t cols, const scalarField& sF)
:
    nRows_(rows),
    nCols_(cols),
    values_(sF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Tensor::~Tensor()
{}


// * * * * * * * * * * * * * * Operator Functions  * * * * * * * * * * * * * //

AFC::scalar AFC::Tensor::operator()(size_t rows, size_t cols) const
{
    return values_[rows * nCols_ + cols];
}


AFC::scalar AFC::Tensor::operator()(size_t rows) const
{
    return operator()(rows, 0);
}


void AFC::Tensor::operator()(size_t rows, size_t cols, scalar value)
{
    values_[rows * nCols_ + cols] = value;
}


void AFC::Tensor::operator()() const
{
    Info<< "\n";
    for (size_t i = 0; i < nRows_; i++)
    {
        Info<< "| ";
        for (size_t j = 0; j < nCols_; j++)
        {
            Info.precision(4);
            Info<< std::scientific << std::setw(13)
                << this->operator()(i, j);
        }
        Info << " |\n";
    }
    Info<< "\n";
}


AFC::scalar& AFC::Tensor::operator()(const size_t rows, const size_t cols)
{
    return values_[rows * nCols_ + cols];
}


AFC::scalar& AFC::Tensor::operator()(const size_t rows)
{
    return operator()(rows, 0);
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


void AFC::Tensor::reset()
{
    for (size_t i = 0; i < nRows_; i++)
    {
        for (size_t j = 0; j < nCols_; j++)
        {
            this->operator()(i, j, scalar(0));
        }
    }
}


// * * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * //


// ************************************************************************* //
