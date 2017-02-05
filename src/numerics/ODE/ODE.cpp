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

AFC::ODE::ODE
(
    const Chemistry& chem
)
{
    if (debug_)
    {
        Info<< "Constructor ODE \n" << endl;
    }

    solver_ = new Seulex(chem); 
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::ODE::~ODE()
{
    if (debug_)
    {
        Info<< "Destructor ODE \n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

void AFC::ODE::solve
(
    const scalar T,
    const scalar p,
    map<word, scalar>& c,
    const scalar dt,
    scalar& dtTry
) const
{
    //- TH::chemistrySolver::ode::solve
    //  {
    //     Jump over
    //     ODESolver::solve
    //     {
    //- dt is the actual time step of the chemistry

    //- Integration starts always with t=0 till dt == tEnd
    //  t represents the actual integration time
    scalar t{0};
    scalar tEnd = dt;

    //- Initialise the state quantities
    StepStatus step(dtTry);

    //- Now we are in the iteration stage
    for (size_t nIter=0; nIter < maxIterations_; nIter++)
    {
        //- Previous time step that we used
        scalar dtTry0 = step.dtTry_;

        step.reject_ = false;

        //- Check if this step is truncated
        //  If yes, set dtTry to the value that we reach tEnd
        //  ATTENSION !!! TODO keep in mind
        if ((t + step.dtTry_) > tEnd)
        {
            step.lastIter_ = true; 
            step.dtTry_ = tEnd - t;
        }

        //- Now we will solve the problem with seulex algorithm
        solver_->solve(T, p, c, t, step);

        //- Checking if we reached the end time and if yes, return the
        //  actual time step that we used
        if (t >= tEnd)
        {
            //- Set dtTry_ to the previous dt step if we made a calculation
            //  and we are in the last step
            if (nIter > 0 && step.lastIter_)
            {
                step.dtTry_ = dtTry0;
            }

            //- The time we return (dtTry is a reference)
            dtTry = step.dtTry_;

            Info<< " Time reached within " << nIter << " loops\n" << endl;

            return;
        }

        //- After first iteration unset bool
        step.firstIter_ = false;

        //- If the dtTry step was rejected by seulex algorithm
        //  we have to go one step back again and treat the step before
        //  differently
        if (step.reject_)
        {
            step.prevReject_ = true;
        }

        //- Next iteration
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
