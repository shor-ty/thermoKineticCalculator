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
#include "numerics.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Numerics::Numerics()
{
    if (debug_)
    {
        Info<< "Constructor Numerics \n" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Numerics::~Numerics()
{
    if (debug_)
    {
        Info<< "Destructor Numerics \n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

AFC::scalar AFC::Numerics::FDMLapacian2ndOrder
(
    const scalar& phi_l,
    const scalar& phi,
    const scalar& phi_r,
    const scalar& Z_l,
    const scalar& Z,
    const scalar& Z_r
) const
{
    //- Return the solution of the laplacian term
    //  l == left position
    //  r == right position
    //  no subscript is the value in the center
    //  Derivation given in [Holzmann]
    return
    (
        //- Enumerator
        (
            (phi_r - phi) / (Z_r - Z)
          - (phi - phi_l) / (Z - Z_l)
        )
        //- Denominator
      / (scalar(0.5) * (Z_r - Z_l))
    );
}

void AFC::Numerics::solveFlamelet
(
    MixtureFraction& flamelet,
    const scalar& chi,
    const scalar& dt
)
{
    //- Save old fields
    const List<map<word, scalar> >& Yold = flamelet.Y();
    const List<map<word, scalar> >& Xold = flamelet.X();
    const List<map<word, scalar> >& Cold = flamelet.C();
    const scalarField& Told = flamelet.T();

    //- Space discretization
    const scalarField& Zvalues = flamelet.Z();
    const int nZ = flamelet.nZPoints();

    //- Species that need to be solved
    const wordList species_ = flamelet.species();

    //- Solve species equation for mass fraction
    //  Time derivation => 1nd Order Euler
    //  Laplace derivation => 2nd Order Linear for non-equal spaced distances
    //  Source term omega handled explicit (for the first implementation)
    //  Equation is solved explicit (implicit should be implemented)

    //- Diffusion coefficient (chi / 2)
    const scalar D = chi/scalar(2);

    //- Source term omega
    const scalar omega = 0;
    
    //- Solve the flamelet equation for each species
    //  Derivation is given in [Holzmann]
    forAll(species_, s)
    {
        //- Temporar field that stores the new values
        scalarField Ynew(nZ, 0);

        //- Copy boundary conditions
        Ynew[0] = Yold[0].at(s);
        Ynew[nZ-1] = Yold[nZ-1].at(s);

        //- Solve for each discrete point (first and last points are BC)
        for (int i = 1; i < nZ-1; ++i)
        {
            Ynew[i] =
                (
                    D
                  * FDMLapacian2ndOrder
                    (
                        Yold[i-1].at(s),
                        Yold[i].at(s),
                        Yold[i+1].at(s),
                        Zvalues[i-1],
                        Zvalues[i],
                        Zvalues[i+1]
                    )
                  + omega
                ) * dt + Yold[i].at(s);
        }

        //- Update the fields (could be done with references faster)
        flamelet.updateY(s, Ynew);
    }

    //- Solve the flamelet equation for the temperature
    {
        //- Temporar field that stores the new temperatures 
        scalarField Tnew(nZ, 0);

        //- Copy boundary conditions
        Tnew[0] = Told[0];
        Tnew[nZ-1] = Told[nZ-1];

        //- Further source terms

        //- Solve the temperature equation
        for (int i = 1; i < nZ-1; ++i)
        {
            Tnew[i] =
                (
                    D
                  * FDMLapacian2ndOrder
                    (
                        Told[i-1],
                        Told[i],
                        Told[i+1],
                        Zvalues[i-1],
                        Zvalues[i],
                        Zvalues[i+1]
                    )
                  + omega
                ) * dt + Told[i];
        }

        //- Update the field
        flamelet.updateT(Tnew);
    }
}


/*
void AFC::Numerics::jacobian
(
    MixtureFraction& mf
)
{
    if (debug)
    {
        Info<< " --> AFC::Numerics::jacobian" << endl;
    }
}
*/

/*void calculate
(
    lookUpTable& lut,
    const scalar& sDR,
    const scalar& defect,
    const unsigned int& nDisPoints
)
{
    //- TODO
    //- Calculate reaction rate factors k for each reaction

    for (unsigned int point=0; point <= nDisPoints; point++)
    {
//        Info<< "Z = " << point << "\n";
        //- Object of discrete mixture fraction
        MixtureFraction& dMF = lut[defect][sDR][point];

        //- Temperature at discrete point Z
        const scalar& T = dMF.T();

        //- Mol fractions at discrete point
//        const map<word, scalar>& speciesMol = dMF.mol();
        {
            //- a) calculate mean molecular weight MW (using mol)
            //dMF.calculateMeanMW("mol");

            //- b) calculate mean heat capacity cp [J/mol/K]
            //dMF.calculateMeanCp(T);

            //- c) calculate mean enthalpy H [J/mol]
            //dMF.calculateMeanH(T);

            //- d) calculate mean entropy S [J/mol/K]
            //dMF.calculateMeanS(T);

            //- e) calculate mean free gibbs energy
            {
             //   const scalar& H = dMF.H();

             //   const scalar& S = dMF.S();

             //   dMF.calculateMeanG(H, S, T);
            }
        }
    }
}

        //- Calculate the source term of species omega
        {
            //- Actual concentration of species at point Zi
            //const map<word, scalar>& con1 = dMF.con();

            //- Copy of species concentration
            //map<word, scalar> con = con1;

            //dMF.calculateOmega(T, con);

            //- Update kf and kb using the new T field

            //- Calc for each species the source term omega
            //chem_.omega(speciesMol);
            //
            //I
        }
    }
}*/





// ************************************************************************* //
