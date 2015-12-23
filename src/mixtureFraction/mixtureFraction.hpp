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

Class
    AFC::MixtureFraction
    
Description
    Abstract AFC::MixtureFraction class for mixtureFraction data and calculation

SourceFiles
    mixtureFraction.cpp

\*---------------------------------------------------------------------------*/

#ifndef MixtureFraction_hpp
#define MixtureFraction_hpp

#include "chemistry.hpp"
#include "thermo.hpp"
#include "transport.hpp"
#include "chemistry.hpp"
#include "properties.hpp"
#include "typedef.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class MixtureFraction Declaration
\*---------------------------------------------------------------------------*/

class MixtureFraction
{
    private:

        // Private Data

            //- Mole fractions X of species at discrete point Z
            map<word, scalar> speciesMol_;

            //- Mass fractions Y of species at discrete point Z
            map<word, scalar> speciesMass_;

            //- Concentration [X] of species at discrete point Z
            map<word, scalar> speciesCon_;

            //- Temperature at discrete point Z
            scalar temperature_{0};

            //- Mean density at discrete point Z
            scalar rho_{0};

            //- Mean molecular weight at discrete point Z
            scalar MW_{0};
        
            //- Density of species i at discrete point Z
            map<word, scalar> rhoSpecies_;

            //- Viscosity at discrete point Z
            scalar mu_{0};

            //- Viscosity of species i at discrete point Z
            map<word, scalar> muSpecies_;

            //- Defect value
            scalar defect_{0};

            //- Mixture fraction value Z
            scalar Z_{0};

    
        // Class data

            //- Thermodynamic class
            const Thermo& thermo_;

            //- Transport class
            const Transport& transport_;

            //- Chemistry class
            const Chemistry& chemistry_;

            //- Properties class
            const Properties& properties_;


        // Debug
        const bool debug{true};

        const bool debugMW{false};


    public:

        //- Constructor
        MixtureFraction
        (
            const Chemistry&,
            const Thermo&,
            const Transport&,
            const Properties&,
            const scalar&,
            const scalar&
        );

        //- Destructor
        ~MixtureFraction();


        // Member functions
        
            //- Calculate the mean molecular weight
            void calcMeanMW
            (
                const word&
            );

            //- Calculate mol fraction out of mass fraction
            void YtoX();

            //- Calculate mass fraction out of mol fraction
            void XtoY();

            //- Calculate concentration out of mass fraction
            void YtoC();

            //- Calculate concentration out of mol fraction
            void XtoC();

            //- Calculate mean density out of mol fraction
            void XtoRho();

            //- Calculate mean density out of mass fraction
            void YtoRho();

            //- Calculate mean density out of concentration fraction
            void CtoRho();
        

        // Check functions


        // Return functions

            void mols(const word) const;

            //- Return species mol fraction
/*            scalar mol
            (
                const word&
            );*/

            //- Return species mol fraction (map)
            map<word, scalar> mol() const;

            //- Tmperature [K]
            scalar T() const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // MixtureFraction_hpp included

// ************************************************************************* //
