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

AFC::Numerics::Numerics
(
//    const int nZ 
)
{
    if (debug)
    {
        Info<< "Constructor Numerics \n" << endl;
    }

    //- Construct matrix for central differencing 
 /*  Info<< "Build matrix\n";

    Info<< "nZ: " << nZ_ << endl;

    M_.resize(nZ_, scalarField(nZ_));

        //- Object of discrete mixture fraction at Z+1 (r: right)
        const MixtureFraction& dMFr = lut[defect][sDR][point+1];

    for (int i=0; i<nZ_; i++)
    {
        for (int j=0; j<nZ_; j++)
        {

            if (i==j)
            {
                M_[i][j] = 2;        
            }
            else if (i==j-1 || i==j+1)
            {
                M_[i][j] = -1;
            }
            else
            {
                M_[i][j] = 0;
            }
        }
    }

    forAll(M_, i)
    {
        forAll(M_[i], j)
        {
            Info<< "  " << M_[i][j];
        }
        Info<< "\n";
    }
    */
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Numerics::~Numerics()
{
    if (debug)
    {
        Info<< "Destructor Numerics \n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //


void AFC::Numerics::solve
(
    scalarField& phi,
    const scalar& dt
)
{
    //- Solve the guy
    /*const int& i = M_.size();

    //-Copy phi
    scalarField phiO = phi;

    for (int Z = 1; Z < i; Z++)
    {
        phi[Z] = phiO[Z] + (phiO[Z-1] - 2*phiO[Z] + phiO[Z+1])/pow(0.05, 2) * dt;
    }

    forAll(phi, Z)
    {
        Info<< Z*0.05 << "  " << phi[Z] << endl; 


                //- c) Update mass fraction of species s using central difference
                {
                    //- Omega in [g/m^3]
                    mass[species[s]] =
                        (
                            sDR/2
                            *(
                                massl.at(species[s])
                              - 2* mass.at(species[s])
                              + massr.at(species[s])
                             )/ pow(dx, 2)
                             + omega.at(species[s]) / dMF.rho()
                        ) * dt + massOld.at(species[s]);
                }
            }
//                std::terminate();
            
            //- d) Update concentration field
            dMF.updateC();*/
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
