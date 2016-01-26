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
            //  [mol/m^3]
            map<word, scalar> speciesCon_;

            //- Temperature at discrete point Z
            scalar temperature_{0};

            //- Mean density at discrete point Z [g/m^3]
            scalar rho_{0};

            //- Mean molecular weight at discrete point Z [g/mol]
            scalar MW_{0};

            //- Mean heat capacity at discrete point Z [J/
            scalar cp_{0};

            //- Mean enthalpy H at discrete point Z
            scalar H_{0};

            //- Mean entropy S at discrete point Z
            scalar S_{0};

            //- Mean free Gibbs energy at discrete point Z
            scalar G_{0};
        
            //- Density of species i at discrete point Z [g/m^3]
            //map<word, scalar> rhoSpecies_;

            //- Viscosity at discrete point Z [m^2/s^2]
            scalar mu_{0};

            //- Viscosity of species i at discrete point Z [m^2/s^2]
            //map<word, scalar> muSpecies_;

            //- Defect value
            scalar defect_{0};

            //- Mixture fraction value Z
            scalar Z_{0};


        // Booleans
        
            //- Mean molecular weight updated?
            bool updatedMW_;


        // Class data
        
            //- Thermodynamic class
            const Thermo& thermo_;

            //- Transport class
            const Transport& transport_;

            // Chemistry class
            Chemistry& chemistry_;

            // Properties class
            const Properties& properties_;


        // Debug
        const bool debug{false};

        const bool debugMW{false};


    public:

        //- Constructor
        MixtureFraction
        (
            Chemistry&,
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
            void calculateMeanMW
            (
                const word& 
            );

            //- Calculate the mean heat capacity
            void calculateMeanCp
            (
                const scalar&
            );

            //- Calculate the mean enthalpy
            void calculateMeanH
            (
                const scalar&
            );

            //- Calculate the mean entropy
            void calculateMeanS
            (
                const scalar&
            );

            //- Calculate the mean free gibs
            void calculateMeanG
            (
                const scalar&
            );

            //- Calculate the mean free gibs
            void calculateMeanG
            (
                const scalar&,
                const scalar&,
                const scalar&
            );

            //- Update the source term omega
            void calculateOmega
            (
                const scalar&,
                map<word, scalar>&
            );

            //- Calculate the formation enthalpy [J/mol]
            void calculateHf
            (
                const word&,
                const scalar&
            ) const;

        // Conversation functions

            //- Calculate mol fraction out of mass fraction
            void YtoX();

            //- Calculate mass fraction out of mol fraction
            void XtoY();

            //- Calculate concentration out of mass fraction
            void YtoC();

            //- Calculate concentration out of mol fraction
            void XtoC();

            //- Calculate mean density out of mol fraction
            void rhoX();

            //- Calculate mean density out of mass fraction
            void rhoY();

            //- Calculate mean density out of concentration fraction
            void rhoC();
        

        // Check functions


        // Return functions

            //- Return species mol fraction (map)
            map<word, scalar> mol() const;

            //- Return species concentration
            map<word, scalar> con();

            //- Return temperature [K]
            scalar T() const;

            //- Return heat capacity [J/mol/K]
            scalar cp() const;

            //- Return heat capacity of species s [J/mol/K]
            scalar calculateCp
            (
                const word&,
                const scalar&
            ) const;

            //- Return mean enthalpy [J/mol]
            scalar H() const;

            //- Return enthalpy of species s [J/mol]
            scalar calculateH
            (
                const word&,
                const scalar&
            ) const;

            //- Return mean entropy [J/mol/K]
            scalar S() const;

            //- Return entropy of species s [J/mol/K]
            scalar calculateS
            (
                const word&,
                const scalar&
            ) const;

            //- Return mean free Gibbs energy
            scalar G() const;

            //- Return free Gibbs energy of species s [J/mol]
            scalar calculateG
            (
                const word&,
                const scalar&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // MixtureFraction_hpp included

// ************************************************************************* //
