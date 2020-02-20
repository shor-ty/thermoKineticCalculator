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
#include "ODE.hpp"
#include "stepStatus.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<typename Type>
AFC::ODE<Type>::ODE(const Chemistry& chem)
:
    StepStatus(1),

    chem_(chem)
{
    solver_ = new Type(chem.species().size());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<typename Type>
AFC::ODE<Type>::~ODE()
{}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

template<typename Type>
void AFC::ODE<Type>::derivative
(
    const scalar dt,
    const scalar T,
    const map<word, scalar>& c0
)
{
    //- Calculate all source terms dc/dt based on
    //  Arrhenius with modificaion of LOW / TROE / SRI / Enhanced
    dcdt_ = chem_.calculateOmega(T, c0);
}


template<typename Type>
void AFC::ODE<Type>::solve
(
    const scalar T,
    const scalar p,
    map<word, scalar>& c,
    const scalar dt,
    scalar& dtTry
)
{
    size_t iterChem{0};

    scalar timeLeft = dt;

    //- Copy of old concentrations
    const map<word, scalar> c0 = c;

    //- TODO SMALL
    while (timeLeft > 1e-15)
    {
        //- Time step of actual chemistry iteration
        scalar dt = 0.00001;

        scalar err = scalar(0);

        //- Calculate the derivative dcdt
        derivative(dt, T, c0);

        //- Loop and adjust ste-size to achieve the desired error
        do
        {
            //- Integration starts always from zero
            solver_->solve(0, c0, dcdt_, dt, c_);

            //- Error is larger than one, reduce dt
            if (err > 1)
            {
                const scalar scale =
                    max(safeScale_*pow(err, -alphaDec_), minScale_);

                dt *= scale;

                //- TODO small
                if (dt < 1e-15)
                {
                    FatalError
                    (
                        "Time step is going to zero",
                        __FILE__,
                        __LINE__
                    );
                }
            }

        }
        while(err > 1);

        timeLeft -= dt;

        iterChem++;
    }
}


template<typename Type>
AFC::scalar AFC::ODE<Type>::normalizeError
(
    const map<word, scalar>& c0,
    const map<word, scalar>& c,
    const map<word, scalar>& err
)
{
    scalar maxErr{0};

    //- TODO
    forAll(err, i)
    {
        //scalar tol = absTol_[i] + relTol_[i] * max(mag(c0[i]), mag(c[i]));
        //maxErr = max(maxErr, mag(err[i])/tol);
        std::terminate();
    }

    return maxErr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
