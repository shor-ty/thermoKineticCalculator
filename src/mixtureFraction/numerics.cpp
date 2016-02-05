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
#include "mixtureFraction.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

// * * * * * * * * * * * * * AFC Numerics functions  * * * * * * * * * * * * //

void calculate
(
    lookUpTable& lut,
    const scalar& sDR,
    const scalar& defect,
    const unsigned int& nPoints,
    const scalar& dt
)
{
    //- Copy of lut
    const lookUpTable lut_old = lut;

    //- Solve the Flamelet-Equation
    //  + Point 0 -> Oxidizer boundary
    //  + Point nPoints -> Fuel boundary
    for (unsigned int point=1; point <= nPoints-1; point++)
    {
        Info<< "Z = " << point << "\n";

        //- Object of discrete mixture fraction
        MixtureFraction& dMF = lut[defect][sDR][point];

        //- Object of old discrete mixture fraction
        const MixtureFraction& dMF_old = lut_old[defect][sDR][point];

        //- Object of discrete mixture fraction at Z-1 (l: left)
        const MixtureFraction& dMFl = lut[defect][sDR][point-1];

        //- Object of discrete mixture fraction at Z+1 (r: right)
        const MixtureFraction& dMFr = lut[defect][sDR][point+1];

        //- Temperature at discrete point Z
        scalar& T = dMF.T();

        //- Species mass fraction calculation
        {
            //- Get all species
            const wordList& species = dMF.species();

            //- Concentration of species at point Z [g/mol]
            map<word, scalar>& con = dMF.con();

            //- Mass fraction of species at point Z 
            map<word, scalar>& mass = dMF.mass();

            //- Mass fraction of species at point Z (old)
            const map<word, scalar>& massOld = dMF_old.mass();

            forAll(species, s)
            {
                Info<< "Concentration of " << species[s] << ": " << con.at(species[s]) << "\n";
            }

            //- Mass fraction of species at Point Z-1
            const map<word, scalar>& massl = dMFl.mass(); 

            //- Mass fraction of species at Point Z+1
            const map<word, scalar>& massr = dMFr.mass();

            //- dx
            const scalar& dx = dMF.Z() - dMFl.Z(); 

            //- Source term of species s
            scalar omega{0};

            //- a) Calculate the source term of species s
            forAll(species, s)
            {
                Info << "++ Update " << species[s] << "\n";
                //- Calculate source term [g/m^3/s]
                omega = dMF.calculateOmega(species[s], T, con);

                //- b) Update density
                {
                    dMF.updateRho(species);
                }

                //- c) Update mass fraction of species s using central difference
                {
                    //- Omega in [g/m^3]
                    Info<< "    " << massl.at(species[s]) << " - " << 2*mass.at(species[s])
                        << " + " << massr.at(species[s]) << " = ";
                    //mass[species[s]] =
                    scalar test =
                        (
                            sDR/2
                            *(
                                massl.at(species[s])
                              - 2* mass.at(species[s])
                              + massr.at(species[s])
                             )/ pow(dx, 2)
                             + omega / dMF.rho()
                        ) * 1e-8 + massOld.at(species[s]);
                    Info<< test << " | " << massOld.at(species[s]) << endl;
                    Info<< "Rho: " << dMF.rho() << "\n";
                    Info<< "omega/rho: " << omega /dMF.rho() << endl;
                }
            }
        
            std::terminate();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
