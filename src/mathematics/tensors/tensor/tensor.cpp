/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     T hermo       | Open Source Thermo-Kinetic Library
  c-o-o-c     K iknetic     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Constructor.

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

\*---------------------------------------------------------------------------*/

#include "definitions.hpp"
#include "tensor.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::Tensor::Tensor()
:
    nRows_(0),
    nCols_(0)
{}


TKC::Tensor::Tensor(const size_t rows, const size_t cols, const scalar value)
:
    nRows_(rows),
    nCols_(cols),
    values_(rows * cols, value)
{}


TKC::Tensor::Tensor(const size_t rows, const size_t cols, const scalarField& sF)
:
    nRows_(rows),
    nCols_(cols),
    values_(sF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::Tensor::~Tensor()
{}


// * * * * * * * * * * * * * * Operator Functions  * * * * * * * * * * * * * //

TKC::scalar TKC::Tensor::operator()(size_t rows, size_t cols) const
{
    return values_[rows * nCols_ + cols];
}


TKC::scalar TKC::Tensor::operator()(size_t rows) const
{
    return operator()(rows, 0);
}


void TKC::Tensor::operator()(size_t rows, size_t cols, scalar value)
{
    values_[rows * nCols_ + cols] = value;
}


void TKC::Tensor::operator()() const
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


TKC::scalar& TKC::Tensor::operator()(const size_t rows, const size_t cols)
{
    return values_[rows * nCols_ + cols];
}


TKC::scalar& TKC::Tensor::operator()(const size_t rows)
{
    return operator()(rows, 0);
}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

size_t TKC::Tensor::cols() const
{
    return nCols_;
}


size_t TKC::Tensor::rows() const
{
    return nRows_;
}


TKC::scalarField TKC::Tensor::values() const
{
    return values_;
}


void TKC::Tensor::reset()
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
