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
            //  [-]
            List<map<word, scalar> > speciesMol_;

            //- Mass fractions Y of species at discrete point Z
            //  [-]
            List<map<word, scalar> > speciesMass_;

            //- Concentration [X] (C) of species at discrete point Z
            //  [mol/m^3]
            List<map<word, scalar> > speciesCon_;

            //- Source term for each species at discrete point Z
            List<map<word, scalar> > omega_;

            //- Concentration [X] of the mixture at discrete point Z
            scalarField Cmix_;

            //- Mixture fraction value Z [-]
            scalarField Z_;

            //- Temperature field for all discrete points Z
            //  [K]
            scalarField T_;

            //- Mean density at discrete point Z
            //  [g/m^3]
            scalarField rho_{0};

            //- Mean molecular weight at discrete point Z
            //  [g/mol]
            scalarField MMW_{0};

            //- Mean heat capacity at discrete point Z
            //  [J/kg/K]
            scalarField cp_{0};

            //- Mean enthalpy H at discrete point Z
            scalarField H_{0};

            //- Mean entropy S at discrete point Z
            scalarField S_{0};

            //- Mean free Gibbs energy at discrete point Z
            scalarField G_{0};

            //- Density of species i at discrete point Z [g/m^3]
            //map<word, scalar> rhoSpecies_;

            //- Thermal conductivity [W/m/K]
            scalarField lambda_{0};

            //- Viscosity at discrete point Z [m^2/s^2]
            scalarField mu_{0};

            //- Viscosity of species i at discrete point Z [m^2/s^2]
            //map<word, scalar> muSpecies_;

            //- Defect value [J/kg]
            scalar defect_{0};

            //- Scalar dissipation rate [1/s]
            scalar sDR_{0};


        // Booleans

            //- Mean molecular weight updated?
            bool updatedMeanMW_{false};

            //- Field that is actual (based on that we choose how to evaluate
            //  the mean molecular weight
            word actualField_;

            //- Mean density updated?
            bool updatedMeanRho_{false};


        // Class data

            //- Thermodynamic class
            const Thermo& thermo_;

            //- Transport class
            const Transport& transport_;

            // Chemistry class
            Chemistry& chemistry_;

            // Properties class
            const Properties& properties_;


    public:

        //- Constructor
        MixtureFraction
        (
            Chemistry&,
            const Thermo&,
            const Transport&,
            const Properties&,
            const scalar,
            const scalar
        );

        //- Destructor
        ~MixtureFraction();


        // Member functions

            //- Calculate density of oxidizer stream
            void calcRhoOxidizer(const Properties&);

            //- Update the mean molecular weight
            void updateMeanMW();

            //- Calculate the mean heat capacity
            void calculateMeanCp(const scalar);

            //- Calculate the mean enthalpy
            void calculateMeanH(const scalar);

            //- Calculate the mean entropy
            void calculateMeanS(const scalar);

            //- Calculate the mean free gibs
            void calculateMeanG(const scalar);

            //- Calculate the mean free gibs
            void calculateMeanG(const scalar, const scalar, const scalar);

            //- Calculate the source term omega for species s
            scalar calculateOmega(const word, const int) const;

            //- Calculate the formation enthalpy [J/mol]
            void calculateHf(const word, const scalar) const;

            //- Output the conservation for the mixture fraction space Z
            void conservation() const;

            //- Output the conservation at discrete point Z
            void conservation(const int) const;

            //- Calculate and return the mass fraction conservation for
            //  the mixture fraction space Z
            scalarField Yconservation() const;

            //- Calculate and return the mass fraction conservation for
            //  the discrete mixture fraction point Z
            scalar Yconservation(const int) const;


        // Update functions

            //- Update the mass fraction Y [-]
            void updateY(const word, const scalarField&);

            //- Update the mass fraction Y [-]
            void updateY(const map<word, scalar>&, const int);

            //- Update the mol fraction x [-]
            void updateX(const map<word, scalar>&, const int);

            //- Update the concentration [x]
            void updateC(const map<word, scalar>&, const int);

            //- Update the mass fraction X [-]
            void updateY();

            //- Update the mol fraction X [-]
            void updateX();

            //- Update the concentration [X] [mol/m^3]
            void updateC();

            //- Update all fields based on the actual field
            void updateFields();

            //- Update the reaction rate coefficients k
            void updatek();

            //- Update the temperature field [K]
            void updateT(const scalarField&);

            //- Update the mean density [kg/m^3]
            void updateRho(const word method = "mass");

            //- Update heat capacity Cp
            void updateCp();


        // Bool switch

            //- Updated mean molecular weight?
            bool updatedMeanMW() const;

            //- Set mean molecular weight switch to false
            void updatedMeanMW(const bool);

            //- Updated mean density
            bool updatedMeanRho() const;

            //- Set mean density switch to false
            void updatedMeanRho(const bool);

            //- Set all bools to false
            void updatedFields(const bool);

            //- Update the source term omega in the mixture fraction space Z
            //  for each species s
            void updateOmega();

            //- Update the source term omega on the discrete mixture fraction
            //  space Z for all species
            void updateOmega(const int);

            //- Update the source term omega on the discrete mixture fraction
            //  space Z and for the species s
            void updateOmega(const int, const word);

            //- Update the source term omega for speacies s
            void updateOmega(const word);



        // Conversation functions

            //- Calculate mol fraction out of mass fraction
            void YtoX();

            //- Calculate mass fraction out of mol fraction
            void XtoY();

            //- Calculate concentration out of mass fraction
            void YtoC();

            //- Calculate concentration out of mol fraction
            void XtoC();

            //- Calculate mass fraction out of concentration
            void CtoY();

            //- Calculate mol fraction out of concentration
            void CtoX();

            //- Calculate mean density out of mol fraction
            void rhoX();

            //- Calculate mean density out of mass fraction
            void rhoY();

            //- Calculate mean density out of concentration fraction
            void rhoC();


        // Return functions

            //- Return the species included in the combustion
            wordList species() const;

            //- Return the species of the oxidizer stream (Z = 0)
            wordList speciesOxidizer() const;

            //- Return the species of the fuel stream (Z = 1)
            wordList speciesFuel() const;

            //- Return number of discrete points
            int nZPoints() const;

            //- Return the value of the mixture fraction Z at position Z
            scalar Z(const int) const;

            //- Return the temperature field at discrete point Z [K] (const)
            scalar T(const int) const;

            //- Return species mol fraction X at discrete point Z (const)
            map<word, scalar> X(const int) const;

            //- Return species mass fraction Y at discrete point Z (const)
            map<word, scalar> Y(const int) const;

            //- Return species concentration [X] at discrete point Z (const)
            map<word, scalar> C(const int) const;

            //- Return species source terms at discrete point Z (const)
            map<word, scalar> omega(const int) const;

            //- Return the scalar dissipation rate of the flamelet [Hz]
            scalar sDR() const;

            //- Return the values of the mixture fraction Z
            scalarField Z() const;

            //- Return the density field [kg/m^3] (const)
            scalarField rho() const;

            //- Return the temperature field [K] (const)
            scalarField T() const;

            //- Return the pressure field [Pa] (const)
            scalar p() const;

            //- Return species mol fraction X for all discrete points Z
            List<map<word, scalar> > X() const;

            //- Return species mass fraction Y for all discrete points Z
            List<map<word, scalar> > Y() const;

            //- Return species concentration [X] for all discrete points Z
            List<map<word, scalar> > C() const;

            //- Return source terms of species for all discrete points Z
            List<map<word, scalar> > omega() const;

            //- Return the mol fraction X for species s for all discr. points
            scalarField X(const word) const;

            //- Return the mass fraction Y for species s for all discr. points
            scalarField Y(const word) const;

            //- Return the concentration [X] for species s for all discr. points
            scalarField C(const word) const;

            //- Return the source term for species s for all discrete points Z
            scalarField omega(const word) const;

            //- Return the mass fraction Y as a map<word, scalarField>
            map<word, scalarField> YField() const;

            //- Return the concentration [X] as a map<word, scalarField>
            map<word, scalarField> CField() const;

            //- Return the concentration [X] for discrete point Z
            map<word, scalar> CField(const int) const;

            //- Return the omega source term field
            map<word, scalarField> omegaField() const;

            //- Return the omega source term field for discrete point Z
            map<word, scalar> omegaField(const int) const;

            //- Return the reaction numbers where species s is included
            List<int> reacNumbers(const word) const;

            //- Return the list of species that are in the elementar reaction r
            wordList speciesInReaction(const int) const;

            //- Return the product species list of reaction r
            wordList speciesProducts(const int) const;

            //- Return the educt species list of reaction r
            wordList speciesEducts(const int) const;

            //- Return the stochiometric factors of the products in reaction r
            map<word, int> nuProducts(const int) const;

            //- Return the stochiometric factors of the educts in reaction r
            map<word, int> nuEducts(const int) const;

            //- Return the elementar reaction r as string
            string elementarReaction (const int) const;

            //- Return the chemistry object
            Chemistry chemistry() const;

            //- Return the thermo object
            Thermo thermo() const;

            //- Return the reaction rate coeff kf of reaction r
            scalar kf(const int, const int) const;

            //- Return the reaction rate coeff kb of reaction r
            scalar kb(const int, const int) const;

            //- Return molecular weight of species s
            scalar MW(const word) const;



        // Write output

            //- Save flamelet tables if required
            void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // MixtureFraction_hpp included

// ************************************************************************* //
