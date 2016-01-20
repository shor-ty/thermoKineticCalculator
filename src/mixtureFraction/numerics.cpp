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
#include "chemistry.hpp"
#include "thermo.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

// * * * * * * * * * * * * * AFC Numerics functions  * * * * * * * * * * * * //


void calculate
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
        //- Object of discrete mixture fraction
        MixtureFraction& dMF = lut[defect][sDR][point];

        //- Temperature at discrete point Z
        const scalar& T = dMF.T();

        //- Mol fractions at discrete point
//        const map<word, scalar>& speciesMol = dMF.mol();

        {
            //- a) calculate mean molecular weight MW (using mol)
            dMF.calculateMeanMW("mol");

            for(int i=1; i<11; i++)
            {
                scalar T = 300 + i*100;

                //- b) calculate mean heat capacity cp [J/mol/K]
                dMF.calculateMeanCp(T);

                //- c) calculate mean enthalpy H [J/mol]
                dMF.calculateMeanH(T);

                //- d) calculate mean entropy S [J/mol/K]
                dMF.calculateMeanS(T);
                
                //- e) calculate mean free gibbs energy
                {
                    const scalar& H = dMF.H();

                    const scalar& S = dMF.S();

                    dMF.calculateMeanG(H, S, T);

                    Info<< dMF.G() << endl;
                }
            }

            //- Calc k with mol fractions and k 
            dMF.k();

            //- Calc for each species the source term omega
            //chem_.omega(speciesMol);
            //
            std::terminate();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
