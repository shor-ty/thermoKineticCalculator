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

#include "transportCalc.hpp"
#include "constants.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::TransportCalc::TransportCalc()
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::TransportCalc::~TransportCalc()
{}


// * * * * * * * * * * * * Reduced collision integrals * * * * * * * * * * * //

AFC::scalar AFC::TransportCalc::reducedCollisionIntegralOmega22
(
    const scalar& Ts
)
{
    //- Constants
    const scalar A = 1.16145;
    const scalar B = 0.14874;
    const scalar C = 0.52487;
    const scalar D = 0.77320;
    const scalar E = 2.16178;
    const scalar F = 2.43787;

    //- Return reduced collision integral omega 2.2
    return (A * pow(Ts, -1*B) + C * (exp(-1*D*Ts)) + E * exp(-1*F*Ts));
}


AFC::scalar AFC::TransportCalc::reducedCollisionIntegralOmega11
(
    const scalar& Ts
)
{
    //- Constants
    const scalar A = 1.06036;
    const scalar B = 0.15610;
    const scalar C = 0.19300;
    const scalar D = 0.47635;
    const scalar E = 1.03587;
    const scalar F = 1.52996;
    const scalar G = 1.76474;
    const scalar H = 3.89411;

    //- Return reduced collision integradl omega 1.1
    return (A / pow(Ts, B) + C / exp(D*Ts) + E / exp(F*Ts) + G / exp(H*Ts));
}


AFC::scalar AFC::TransportCalc::rho
(
    const scalar& MW, 
    const scalar& T,
    const scalar& p 
)
{
    //- Return density [kg/m^3]
    return (p * MW / AFC::Constants::R / T); 
}


// * * * * * * * * Calculation Functions For Viscosity  * * * * * * * * * * *//

void AFC::TransportCalc::viscosity
(
    const Thermo& thermo,
    TransportData& transData
)
{
    Info<< "    c-o Pure Species viscosity calculation ... ";

    //- Get species used in kinetic calculations
    const wordList& chemSpecies = transData.chemSpecies();

    //- Calculate viscosity values of chemical species
    forAll(chemSpecies, species)
    {
        //- Methode of Chung et. al. 1984, 1988

        //- Methode of Neufeld et. al. 1972 [Pa s]
        scalar T = 1000;
        const scalar& visco = viscosityNeufeld(species, T, thermo, transData);

        Info<< species << ": " << visco << " [Pa s] \n";
    }
}


AFC::scalar AFC::TransportCalc::viscosityNeufeld
(
    const word& species,
    const scalar& T,
    const Thermo& thermo,
    const TransportData& transData
)
{
    //- Viscosity calculation suggested by Neufeld et. al. (1972)

    //- Molecular weight [g/mol]
    const scalar& MW = thermo.MW(species);

    //- Lennard-Jones potential well depth eps/kb [K]
    const scalar& LJP = transData.LJP(species);

    //- Lennard-Jones collision diameter [Angstroms]
    const scalar& sigma = transData.LJCD(species);

    //- Dimensionless temperature T*
    const scalar Ts =   T * pow(LJP, -1);

    //- Valid for 0.3 <= Ts <= 100
    if (Ts < 0.3 || Ts > 100)
    {
        //- Output warning to file
    }

    //- Calculate collision integral omega 22
    const scalar omegaColl = reducedCollisionIntegralOmega22(Ts);

    //- Calculate viscosity in micro Poise (P) 
    //  [P] = 0.1 [Pa s]
    const scalar mu = 26.693 * sqrt(MW * T) / (pow(sigma, 2) * omegaColl);

    //- Return vicsosity as [Pa s] = [kg/m/s]
    return mu * 1e-7;
}


void AFC::TransportCalc::viscosityChung
(
    const word& species
)
{
    // Constants

    //- Stefan Bolzmann constant [J/m^2/K^4]
    /*const scalar& kB = AFC::Constants::stefanBolzmann;

    //- Lennard-Jones collision diameter [Angstroms]
    const scalar& LJCD = transData.LJCD(species);

    //- Lennard-Jones potential well depth eps/kb [K]
    const scalar& LJP = transData.LJP(species);

    //- Dipole moment [debey]
    const scalar& muk = transData.muk(species);

    //- Moleculare weight [g/mol]
    const scalar& MW = thermo.MW(species);

    //- Pi [-]
    const scalar& pi = M_PI;

    //- a) Reduced dipole moment
    const scalar delta = 0.5 * muk * muk / (LJP * pow(LJCD, 3));

    //- Start temperature [K]
    scalar T{300};

    //- c) 
    for(; T < 4000.;)
    {

        //- Dimensionless temperature T*
        const scalar Ts = 1.2593 * T

        const scalar mu =
            5. * sqrt(pi * MW * kB * T)
         / (16. * pi * LJCD * LJCD);


        std::terminate();
        T += 50;
    }*/
}


// * * * * * * * * Calculation Functions For Conducitvity  * * * * * * * * * //

void AFC::TransportCalc::thermalConductivity
(
    const Thermo& thermo,
    TransportData& transData
)
{
    Info<< "    c-o Pure Species thermal conductivity calculation ... ";

    //- Get species used in kinetic calculations
    const wordList& chemSpecies = transData.chemSpecies();

    //- Calculate viscosity values of chemical species
    forAll(chemSpecies, species)
    {
        //- Methode mentioned by Warnatz
        scalar T = 1000;
        const scalar& thermalConductivity = 
            thermalConductivityWarnatz(species, T, thermo, transData);

        Info<< species << ": " << thermalConductivity << " [W/m/K]\n";
    }
}


AFC::scalar AFC::TransportCalc::thermalConductivityWarnatz3
(
    const word& species,
    const scalar& T,
    const Thermo& thermo,
    const TransportData& transData
)
{
    //- Calculation of thermal conductivity by Warnatz

    //- Molecular weight [g/mol]
    /*const scalar& MW = thermo.MW(species);

    //- Universal gas constant [J/mol/K]
    const scalar& R = AFC::Constants::R;

    //- Lennard-Jones potential well depth eps/kb [K]
    const scalar& LJP = transData.LJP(species);

    //- Lennard-Jones collision diameter [Angstroms]
    const scalar& sigma = transData.LJCD(species);

    //- Dimensionless temperature T*
    const scalar Ts =   T * pow(LJP, -1);

    //- Linear molecule
    {
        //- Calc CvTrans
        const scalar CvTrans = 3./2. * R;

        //- Calc CvRot
        const scalar CvRot = 1 * R;

        //- Calc CvVib
        const scalar CvVib = Cv - 5./2. * R
    }

    //- Non-Linear molecule
    {
        //- Calc CvTrans
        const scalar CvTrans = 3./2. * R;

        //- Calc CvRot
        const scalar CvRot = 3./2. * R;

        //- Calc CvVib
        const scalar CvVib = Cv - 3. * R
    }

    //- Calculate A
    const scalar A = 5./2. - rho * Dkk / mu;

    //- Calculate B
    const scalar B = ZRot + 2/M_PI * (5./3. * CvRot / R + rho * Dkk / mu);

    //- Calculate fTrans
    const scalar fTrans = 5./2. * (1 - (2*CvRot * A) / (M_PI * CvTrans * B);

    //- Calculate thermal conductivity 
    const scalar thermalConductivity =
        mu/MW * (fTrans * CvTrans + fRot * CvRot + fVib * CvVib);
    */
    return T;
    
}


AFC::scalar AFC::TransportCalc::thermalConductivityWarnatz
(
    const word& species,
    const scalar& T,
    const Thermo& thermo,
    const TransportData& transData
)
{
    //- Calculation of thermal conductivity by Warnatz

    //- Molecular weight [kg/mol]
    const scalar MW = thermo.MW(species)/1000;

    //- Lennard-Jones potential well depth eps/kb [K]
    const scalar& LJP = transData.LJP(species);

    //- Lennard-Jones collision diameter [Angstroms]
    const scalar& sigma = transData.LJCD(species);

    //- Dimensionless temperature T*
    const scalar Ts =   T * pow(LJP, -1);

    //- Valid for 0.3 <= Ts <= 100
    if (Ts < 0.3 || Ts > 100)
    {
        //- Output warning to file
    }

    //- Calculate collision integral omega 22
    const scalar omegaColl = reducedCollisionIntegralOmega22(Ts);

    //- Calculate thermal conductivity in [J/cm/K/s]
    //const scalar lambda = 
    //    8.323e-6 * sqrt(T/MW) / (pow(sigma, 2) * omegaColl);
    const scalar lambda = 8.323e-6 * sqrt(T / MW) / (pow(sigma * 0.1, 2) * omegaColl);

    //- Return thermal conductivity in [W/m/K]
    return lambda;
}
    

// * * * * * * * Calculation Functions For Binary Diffusivity * * * * * * * *//

void AFC::TransportCalc::binaryDiffusivity
(
    const Thermo& thermo,
    TransportData& transData
)
{
    Info<< "    c-o Pure Species binary diffusivity calculation ... ";

    //- Get species used in kinetic calculations
    const wordList& chemSpecies = transData.chemSpecies();

    //- Generate species combinations
    map<word, wordList> combinations;

    forEach(chemSpecies, s)
    {
        for(unsigned i = s+1; i < chemSpecies.size(); i++)
        {
            Info<< "Species: " << chemSpecies[s] << " ... " << chemSpecies[i] << "\n";
            combinations[chemSpecies[s]].push_back(chemSpecies[i]); 
        }
    }

    //- Pressure [Bar]
    const scalar& p = thermo.p()/1e5;

    Info<< "Pressure: " << p << endl;

    const scalar T = 300;

    //- Loop through the combinations
    forEach(chemSpecies, s)
    {
        //- Last species break
        if (s == chemSpecies.size()-1)
        {
            break;
        }

        //- First species
        const word& species1 = chemSpecies[s];

        const wordList& pairs = combinations.at(species1);

        forEach(pairs, i)
        {
            //- Second species
            const word& species2 = pairs[i];

            //- Molecular weight
            const scalar& MW1 = thermo.MW(species1);
            const scalar& MW2 = thermo.MW(species2);

            //- Lennard-Jones potential well depth eps/kb [K]
            const scalar& LJP1 = transData.LJP(species1);
            const scalar& LJP2 = transData.LJP(species2);

            //- Lennard-Jones collision diameter [Angstroms]
            const scalar& sigma1 = transData.LJCD(species1);
            const scalar& sigma2 = transData.LJCD(species2);

            //- Combined parameters 
            const scalar& MW12 = 2 * pow((1/MW1 + 1/MW2), -1);
            const scalar& LJP12 = sqrt(LJP1 * LJP2);
            const scalar& sigma12 = (sigma1 + sigma2) / 2;

            //- Reduced temperature Ts12
            const scalar& Ts12 = T * pow(LJP12, -1);

            //- Calculate collision integral omega 22
            const scalar omegaColl = reducedCollisionIntegralOmega11(Ts12);

            //- Calculate binary diffusivity [cm^2/s]
            const scalar& Dkk =
                0.002662 * sqrt(pow(T, 3))
              / (p * sqrt(MW12) * pow(sigma12, 2) * omegaColl);

            //Info << sqrt(pow(T, 3)) << "/" << "(1 * " << sqrt(MW12) << " * " << pow(sigma12,2)
             //   << " * "<< omegaColl << ")\n";

            //Info<< (p * sqrt(MW12) * pow(sigma12, 2) * omegaColl) << "\n";

            Info << s+1 << "   " << i+1 << "    " 
                << species1 << "   " << species2 << "    " << Dkk/10000 << endl;
        }
    }
}


// ************************************************************************* //
