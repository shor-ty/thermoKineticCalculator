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
#include "time.hpp"
#include "timeReader.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::Time::Time(const string path)
{
    TimeReader reader(path);

    reader.read(*this);

    checkData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::Time::~Time()
{}


// * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * //

void TKC::Time::operator+=(const scalar deltaT)
{
    Info<< "Operator" << runTime() << " + " << deltaT << endl;
    runTime_ = runTime() + deltaT;
}


// * * * * * * * * * * * * * * * Insert Functions  * * * * * * * * * * * * * //

void TKC::Time::endTime(const scalar value)
{
    endTime_ = value;
}


void TKC::Time::writeTime(const scalar value)
{
    writeTime_ = value;
}


void TKC::Time::dTKinetic(const scalar value)
{
    deltaTKinetic_ = value;
}


void TKC::Time::dTFlow(const scalar value)
{
    deltaTFlow_ = value;
}


// * * * * * * * * * * * * * * * Return function * * * * * * * * * * * * * * //

const TKC::scalar TKC::Time::runTime() const
{
    return runTime_;
}


const TKC::scalar TKC::Time::endTime() const
{
    return endTime_;
}


const TKC::scalar TKC::Time::writeTime() const
{
    return writeTime_;
}


const TKC::scalar TKC::Time::dTKinetic() const
{
    return deltaTKinetic_;
}


const TKC::scalar TKC::Time::dTFlow() const
{
    return deltaTFlow_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const bool TKC::Time::loop() const
{
    if (runTime_ == endTime_)
    {
        return 0;
    }
    else if (runTime_ < endTime_)
    {
        return 1;
    }
    else
    {
        Warning
        (
            "In run time is larger than the end time of the calculation",
            __FILE__,
            __LINE__
        );

        return 0;
    }
}


void TKC::Time::checkData() const
{
    //- Do some checks if input makes sense or something is missing
    if (deltaTKinetic_ == 0)
    {
        ErrorMsg
        (
            "The initial time step for the kinetic calculation is zero",
            __FILE__,
            __LINE__
        );
    }

    if (deltaTFlow_ == 0)
    {
        ErrorMsg
        (
            "The initial time step for the flow calculation is zero",
            __FILE__,
            __LINE__
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
