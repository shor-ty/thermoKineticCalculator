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
#include "seulex.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Seulex::Seulex
(
    const Chemistry& chem
)
:
    //- Rate of concentration change
    //dcdt_("H2", 0),

    //- Square matrix for Jacobian
    dcdc_
    (
        chem.species().size(),
        chem.species().size()
    )

    //- TODO relTol (min(1e-4, relTol_))

//    nSeq_(iMax_),

//    jacRedo_(1e-5),

//    theta_(2*jacRedo_),

//    cpu_(iMax_),

//    coeff_(iMax_, iMax_)

//    table_(kMax_, n_),

//    dfdt_(n_),

//    dfdc_(n_)
{
    if (debug_)
    {
        Info<< "Constructor Seulex \n" << endl;
    }

    //- Make new object for Jacobion calculation
    jac_ = new Jacobian(chem);

    //- Init the dcdt field
    {
        const wordList& species = chem.species();

        forAll(species, s)
        {
            dcdt_[s] = scalar(0);
        }
    }

    /*
    //- Quantities to evaluate the factors for the major parts of the algorithm
    const scalar cpuFunc{1};

    const scalar cpuJac{5};

    const scalar cpuLU{1};

    const scalar cpuSolve{1};

    nSeq_[0] = 2;
    nSeq_[1] = 3;

    for (size_t i=2; i<iMax_; i++)
    {
        nSeq_[i] = 2*nSeq_[i-2];
    }

    cpu_[0] = cpuJac + cpuLU + nSeq_[0]*(cpuFunc + cpuSolve);

    for (size_t k=0; k<kMax_; k++)
    {
        cpu_[k+1] = cpu_[k] + (nSeq_[k+1]-1)*(cpuFunc + cpuSolve) + cpuLU;
    }

    // Calculate and set the extrapolation coeff array
    for (size_t k=0; k<iMax_; k++)
    {
        for (size_t l=0; l<k; l++)
        {
            const scalar ratio = scalar(nSeq_[k])/nSeq_[l];

            coeff_(k, l, 1/(ratio -1));
        }
    }
    */
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Seulex::~Seulex()
{
    if (debug_)
    {
        Info<< "Destructor Seulex \n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

void AFC::Seulex::solve
(
    const scalar T,
    const scalar p,
    map<word, scalar>& c,
    scalar& t,
    StepStatus& step 
) const
{
//    temp_[0] = 1e15;

    scalar dt = step.dtTry_;

/*    c0_ = c;

    dtOpt_[0] = fabs(0.1*dt);

    if (step.firstIter_ || step.prevReject_)
    {
        theta_ = 2*jacRedo_;
    }

    if (step.firstIter_)
    {
        //- TODO relTol
//        const scalar logTol = -log10(10-5 + 1e-8) * 0.6 + 0.5;
//        kTarget_ = max(1, min(kMax_ - 1, int(logTol)));
    }

    bool jacUpdated = false;
    */

    //if (theta_ > jacRedo_)
    //{
        jac_->jacobian(T, p, t, c, dcdt_, dcdc_);
    //}
    //
    forAll(dcdt_, t)
    {
        Info<< t.first << "  " << t.second << endl;
    }

    dcdc_();
    std::terminate();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
