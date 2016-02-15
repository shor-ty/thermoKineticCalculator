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
        //Info<< " Z = " << point << "\n";

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

            //- Mass fraction of species at Point Z-1
            const map<word, scalar>& massl = dMFl.mass(); 
            const scalar& Tl = dMFl.T();
            const scalar& cpl = dMFl.cp();

            //- Mass fraction of species at Point Z+1
            const map<word, scalar>& massr = dMFr.mass();
            const scalar& Tr = dMFr.T();
            const scalar& cpr = dMFr.cp();

            //- dx
            const scalar& dx = dMF.Z() - dMFl.Z(); 

            //- Source term of species s
            map<word, scalar> omega;

            //- a) Calculate the source term of species s
            forAll(species, s)
            {
                Info << "  ++ Update " << species[s] << " --> " << con.at(species[s]) << "\n";

                //- Calculate source term [g/m^3/s]
                omega[species[s]] = dMF.calculateOmega(species[s], T, con);


                //- b) Update density
                {
                    dMF.updateC();
                    dMF.updateRho();
                }

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
            dMF.updateC();

            //- e) Update density
            dMF.updateRho();

            //- f) Update heat capacity
            dMF.updateCp();

            //- g) Calculate source term for temperature
            scalar source{0};

            forAll(species, s)
            {
                source += dMF.calculateH(species[s], T) * omega.at(species[s]);
            }
        
            Info<< "SDR:: " << sDR << endl;
            Info<< "cp:: " << dMF.cp() << endl;
            Info<< "cpl:: " << cpl << endl;
            Info<< "cpr:: " << cpr << endl;
            //- h) Updating Temperature using central difference
            scalar tT = sDR/2 * ( Tl - 2 * T + Tr) / pow(dx, 2) - source
                + sDR / ( 2 * dMF.cp() ) * (( cpr - cpl ) / (2*dx)) ; 
            
            Info<< "Temperature = " << tT << endl;
            


//            std::terminate();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
