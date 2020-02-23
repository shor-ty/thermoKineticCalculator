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
#include "matrix.hpp"
#include "vector.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::TransportCalc::TransportCalc()
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::TransportCalc::~TransportCalc()
{}


// * * * * * * * * * * * * Reduced collision integrals * * * * * * * * * * * //

AFC::scalar
AFC::TransportCalc::reducedCollisionIntegralOmega22(const scalar Ts) const
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


AFC::scalar
AFC::TransportCalc::reducedCollisionIntegralOmega11(const scalar Ts) const
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


AFC::scalar
AFC::TransportCalc::rho(const scalar MW, const scalar T, const scalar p) const
{
    //- Return density [kg/m^3]
    return (p * MW / AFC::Constants::R / T); 
}


// * * * * * * * * Calculation Functions For Viscosity  * * * * * * * * * * *//

AFC::scalar AFC::TransportCalc::viscosity
(
    const word species,
    const scalar T,
    const Thermo& thermo,
    const TransportData& transData,
    const word method
) const
{
    //- Molecular weight [kg/mol]
    const scalar MW = thermo.MW(species);

    //- Lennard-Jones potential well depth eps/kb [K]
    const scalar LJP = transData.LJP(species);

    //- Lennard-Jones collision diameter [Angstroms]
    const scalar sigma = transData.LJCD(species);

    //- Method of Hirschfelder et. al. 1954 [kg/m/s]
    if (method == "Hirschfelder")
    {
        return viscosityHirschfelder(species, T, MW, LJP, sigma);
    }

    //- Method of Chung et. al. 1984, 1988
    else if (method == "Chung")
    {
        NotImplemented
        (
            __FILE__,
            __LINE__
        );

        //- For compiler
        return -1; 
    }

    //- Polynomial fit
    else if (method == "Polynomial")
    {
        const scalarField& polyCoeffs =
            transData.viscosityPolyCoeffs(species); 

        return viscosityPolynomial(T, polyCoeffs);
    }
    else
    {
        ErrorMsg
        (
            "   Method " + method + " is not available for calculation of\n"
            "   viscosity. You have to implement it yourself or ask for help",
            __FILE__,
            __LINE__
        );

        //- For compiler
        return -1;
    }
}


AFC::scalar AFC::TransportCalc::viscosityHirschfelder
(
    const word species,
    const scalar T,
    const scalar MW,
    const scalar LJP,
    const scalar sigma
) const
{
    //- Viscosity calculation suggested by Hirschfelder et. al. (1954)

    //- Dimensionless temperature T*
    const scalar Ts =   T * pow(LJP, -1);

    //- Valid for 0.3 <= Ts <= 100
    if (Ts < 0.3 || Ts > 100)
    {
        //- Output warning to file
    }

    //- Calculate collision integral omega 22
    const scalar omegaColl = reducedCollisionIntegralOmega22(Ts);

    //- Calculate viscosity in micro Poise 10^-6 [P]
    //  [P] = 0.1 [Pa s]
    //  MW [kg/mol] here we need [g/mol]
    const scalar mu = 26.693 * sqrt(MW * 1e3 * T) / (pow(sigma, 2) * omegaColl);

    //- Return vicsosity as [Pa s] = [kg/m/s]
    return (mu * 1e-7);
}


/*AFC::scalar AFC::TransportCalc::viscosityChung
(
    const word& species,
    const scalar& T,
    const Thermo& thermo,
    const TransportData& transData
) const
{
    // Constants

    //- Stefan Boltzmann constant [J/m^2/K^4]
    const scalar& kB = AFC::Constants::kB;

    //- Lennard-Jones collision diameter [Angstroms]
    const scalar& LJCD = transData.LJCD(species);

    //- Lennard-Jones potential well depth eps/kb [K]
    const scalar& LJP = transData.LJP(species);

    //- Dipole moment [debey]
    const scalar& muk = transData.muk(species);

    //- Moleculare weight [g/mol]
    const scalar& MW = thermo.MW(species);

    //- Rotational relaxation collision number Zrot at 298K
    const scalar& Zrot298 = transData.ZRot298(species);

    //- Pi [-]
    const scalar& pi = M_PI;

    //- a) Reduced dipole moment
    const scalar delta = 0.5 * muk * muk / (LJP * pow(LJCD, 3));

    //- Geometric configuration TODO
    const scalar geometricConfig = 1; //= transData.geometricConfig(species);

    //- Mono atomic
    if (geometricConfig == scalar(0))
    {
        //- Lambda [erg/cm/K/s]
        const scalar lambda = 
            viscosity(species, T, thermo, transData)
          / MW * scalar(15)/scalar(4) * AFC::Constants::Rerg;

        //- Return lambda
        return lambda;
    }
    else
    {
        scalar Cvtrans{0};
        scalar Cvrot{0};
        scalar Cvvib{0};

        //- Linear atomic
        if (geometricConfig == scalar(1))
        {
            Cvtrans = scalar(3)/scalar(2) * AFC::Constants::Rerg;

            Cvrot = AFC::Constants::Rerg;

            Cvvib = 
                thermo.cv(species, T)
              - scalar(5)/scalar(2) * AFC::Constants::Rerg;
        }
        else if (geometricConfig == scalar(2))
        {
            Cvtrans = scalar(3)/scalar(2) * AFC::Constants::Rerg;

            Cvrot = Cvtrans;

            Cvvib = thermo.cv(species, T) - 3 * AFC::Constants::Rerg;
        }

        //- TODO units
        const scalar& nu = viscosity(species, T, thermo, transData);
        //const scalar& Dii = binaryDiffusivity(species, T, thermo, transData);

        //- Denisty in [g/m^3]
        //const scalar& rho = thermo.rho(species, T);
        
        //const scalar& A = scalar(5)/scalar(2) - rho * Dii / nu;
        //const scalar& B =
        //    Zrot298 + scalar(2) / pi
        // * (scalar(5)/scalar(3) * Cvrot / AFC::Constants::Rerg + rho * Dii / nu );

        //const scalar ftrans =
        //    scalar(5)/scalar(2) * (1 - 2/pi * Cvrot * A / Cvtrans / B);

        //const scalar frot = rho * Dii / nu * (1 + 2 * A / pi / B);

        //const scalar fvib = rho * Dii / nu;
        return 0;
    }
}*/


AFC::scalar AFC::TransportCalc::viscosityPolynomial
(
    const scalar T,
    const scalarField& polyCoeffs 
) const
{
    return
    (
        exp
        (
            polyCoeffs[3]
          + polyCoeffs[2] * log(T)
          + polyCoeffs[1] * pow(log(T), 2)
          + polyCoeffs[0] * pow(log(T), 3)
        )
    );
}


void AFC::TransportCalc::fitViscosity
(
    TransportData& transData,
    const Thermo& thermo
) const
{
    //- TODO use one function instead of 3
    //- Species from chemistry
    const wordList& species = transData.chemistrySpecies();

    //- Fit for all species
    forAll(species, s)
    {
        //- Temperature field 
        scalarField T;

        //- Mu field
        scalarField mu;

        //- Fitting curve are made out of n discrete points from 300K - 4300K
        //  TODO bounds as parameters
        const unsigned int n{10};

        const scalar dT = scalar(4000)/n;

        //- Fill the fields with necessary data
        for (size_t i = 0; i < n; ++i)
        {
            //- Temperature value
            const scalar Tv = scalar(300) + i * dT;

            //- Add temperature
            T.push_back(Tv);

            //- Add viscosity
            mu.push_back(log(viscosity(s, Tv, thermo, transData)));
        }

        //- Fitting 
        //  This could be done more beatuiful with objects
        //
        //  Polynom: ln(mu) = D ln(T)^3 + C ln(T)^2  + B ln(T) + A 
        //  f1(T) = ln(T)^3
        //  f2(T) = ln(T)^2
        //  f3(T) = ln(T)
        //  f4(T) = 1

        //- Now we need to build a matrix A that looks like 
        //          1             2             3             4          5
        //1 | [f1(T)f1(T)]  [f1(T)f2(T)]  [f1(T)f3(T)]  [f1(T)f4(T)] [yf1(T)] |
        //2 | [f2(T)f1(T)]  [f2(T)f2(T)]  [f2(T)f3(T)]  [f2(T)f4(T)] [yf1(T)] |
        //3 | [f3(T)f1(T)]  [f3(T)f2(T)]  [f3(T)f3(T)]  [f3(T)f4(T)] [yf1(T)] |
        //4 | [f4(T)f1(T)]  [f4(T)f2(T)]  [f4(T)f3(T)]  [f4(T)f4(T)] [yf1(T)] |
        //
        //  [] means gauss summation; y = mu value
        //  Matrix is symmetric in case of the first 4 rows and cols
        scalar gS_11{0};
        scalar gS_12{0};
        scalar gS_13{0};
        scalar gS_14{0};
        scalar gS_15{0};

        scalar gS_22{0};
        scalar gS_23{0};
        scalar gS_24{0};
        scalar gS_25{0};

        scalar gS_33{0};
        scalar gS_34{0};
        scalar gS_35{0};

        scalar gS_44{0};
        scalar gS_45{0};

        //- a) Calculate the gauss summations
        for (size_t j = 0; j < n; ++j)
        {
            // Temperature value
            const scalar& Tv = T[j];
            const scalar& muv = mu[j];

            gS_11 += pow(log(Tv), 3) * pow(log(Tv), 3);
            gS_12 += pow(log(Tv), 3) * pow(log(Tv), 2);
            gS_13 += pow(log(Tv), 3) * pow(log(Tv), 1);
            gS_14 += pow(log(Tv), 3) * 1;
            gS_15 += muv * pow(log(Tv), 3);

            gS_22 += pow(log(Tv), 2) * pow(log(Tv), 2);
            gS_23 += pow(log(Tv), 2) * pow(log(Tv), 1);
            gS_24 += pow(log(Tv), 2) * 1;
            gS_25 += muv * pow(log(Tv), 2);

            gS_33 += pow(log(Tv), 1) * pow(log(Tv), 1);
            gS_34 += pow(log(Tv), 1) * 1;
            gS_35 += muv * pow(log(Tv), 1);

            gS_44 += 1 * 1;
            gS_45 += muv * 1;
        }

        //- b) Assign the values to the matrix A
        Matrix A(4,4);

        A(0,0,gS_11);
        A(0,1,gS_12);
        A(0,2,gS_13);
        A(0,3,gS_14);

        //- A_21 = A_12
        A(1,0,gS_12);
        A(1,1,gS_22);
        A(1,2,gS_23);
        A(1,3,gS_24);

        //  A_31 = A_13
        //  A_32 = A_23
        A(2,0,gS_13);
        A(2,1,gS_23);
        A(2,2,gS_33);
        A(2,3,gS_34);

        //  A_41 = A_14
        //  A_42 = A_24
        //  A_43 = A_34
        A(3,0,gS_14);
        A(3,1,gS_24);
        A(3,2,gS_34);
        A(3,3,gS_44);

        //- Calculate the inverse
        Matrix AI = A.inverse();

        //- Init the solution vector
        Vector b(4);

        b(0, gS_15);
        b(1, gS_25);
        b(2, gS_35);
        b(3, gS_45);

        //- Calculate the polynomial coefficients
        //  x = b * A^-1
        const Vector x = b * AI;

        //- Assign the values
        transData.viscosityPolyCoeffs(s, x);
    }
}


// * * * * * * * * Calculation Functions For Conducitvity  * * * * * * * * * //

AFC::scalar AFC::TransportCalc::thermalConductivity
(
    const word species,
    const scalar T,
    const Thermo& thermo,
    const TransportData& transData,
    const word method
) const
{
    //- Methode mentioned by Warnatz
    if (method == "Warnatz")
    {
        return (thermalConductivityWarnatz(species, T, thermo, transData));
    }
    //- Methode mentioned by Warnatz and given in Chemkin Collection
    if (method == "WarnatzCC")
    {
        NotImplemented
        (
            __FILE__,
            __LINE__
        );

        //return (thermalConductivityWarnatzCC(species, T, thermo, transData));
        return -1;
    }
    else if (method == "Polynomial")
    {
        //- Polynomial coefficients
        const scalarField& polyCoeffs =
            transData.thermalConductivityPolyCoeffs(species);

        return (thermalConductivityPolynomial(T, polyCoeffs));
    }
    else
    {
        ErrorMsg
        (
            "    Method " + method + " is not available for calculation of\n"
            "    thermal conductivity. You have to implement it yourself\n"
            "    or ask for help",
            __FILE__,
            __LINE__
        );

        //- For compiler
        return -1;
    }
}


AFC::scalar AFC::TransportCalc::thermalConductivityWarnatz
(
    const word species,
    const scalar T,
    const Thermo& thermo,
    const TransportData& transData
) const
{
    //- Calculation of thermal conductivity by Warnatz

    //  Molecular weight [kg/mol]
    const scalar MW = thermo.MW(species);

    //- Lennard-Jones potential well depth eps/kb [K]
    const scalar LJP = transData.LJP(species);

    //- Lennard-Jones collision diameter [Angstroms]
    //  We need [m] in the formula -> 1 [Anstroms] = 1e-10 [m]
    const scalar sigma = transData.LJCD(species) * scalar(1e-10);

    //- Dimensionless temperature T*
    const scalar Ts = T * pow(LJP, -1);

    //- Valid for 0.3 <= Ts <= 100
    if (Ts < 0.3 || Ts > 100)
    {
        //- Output warning to file
        Warning
        (
            "      Calculation of thermal conductivity based on Warnatz is\n"
            "      not in the valid range for species " + species + " anymore."
            "\n      Valid range is:  0.3 < Ts < 100 while Ts = "
            + toStr(Ts) + "\n      Going on in calculation...\n",
            __FILE__,
            __LINE__
        );
    }

    //- Calculate collision integral omega 22
    const scalar omegaColl = reducedCollisionIntegralOmega22(Ts);

    //- Calculate thermal conductivity in [J/cm/K/s]
    //const scalar lambda = 
    //    8.323e-6 * sqrt(T/MW) / (pow(sigma, 2) * omegaColl);

    //- Calculate thermal conductivity in [J/cm/K/s]
    //const scalar lambda =
    //    8.323e-6 * sqrt(T / MW) / (pow(sigma * 0.1, 2) * omegaColl);
    
    //- Calculate thermal conductivity in [W/m/K]
    //  [Poling]
    const scalar lambda =
        2.63 * 1e-23 * sqrt(T / MW) / ( pow(sigma, 2) * omegaColl );

    //- Return thermal conductivity in [W/m/K]
    return lambda;
}


/*AFC::scalar AFC::TransportCalc::thermalConductivityWarnatzCC
(
    const word& species,
    const scalar& T,
    const Thermo& thermo,
    const TransportData& transData
) const
{
    //- Calculation of thermal conductivity by Warnatz

    //- Molecular weight [g/mol]
    const scalar& MW = thermo.MW(species);

    //- Universal gas constant [J/mol/K]
    const scalar& R = AFC::Constants::R;

    //- Lennard-Jones potential well depth eps/kb [K]
    const scalar& LJP = transData.LJP(species);

    //- Lennard-Jones collision diameter [Angstroms]
    const scalar& sigma = transData.LJCD(species);

    //- Dimensionless temperature T*
    const scalar Ts =   T * pow(LJP, -1);

    //- Constant heat capacity at constant volume
    const scalar Cv = thermo.cv(species, T);
    
    //- Geometric configuration
    const int geometricConfig = transData.geometricalConfig(species);

    scalar CvTrans{0};
    scalar CvRot{0};
    scalar CvVib{0};

    if (geometricConfig == 0)
    {
        CvTrans = 3./2. * R; 
    } 
    //- Linear molecule
    else if (geometricConfig == 1)
    {
        //- Calc CvTrans
        CvTrans = 3./2. * R;

        //- Calc CvRot
        CvRot = 1 * R;

        //- Calc CvVib
        CvVib = Cv - 5./2. * R;
    }
    //- Non-Linear molecule
    else if (geometricConfig == 2)
    {
        //- Calc CvTrans
        CvTrans = 3./2. * R;

        //- Calc CvRot
        CvRot = 3./2. * R;

        //- Calc CvVib
        CvVib = Cv - 3. * R;
    }

    //- Binary Diffusivity Dkk [m^2/s]
    const scalar Dkk =
        binaryDiffusivity(species, species, T, thermo, transData);
    
    //- Viscosity mu [kg/m/s]
    const scalar mu = viscosity(species, T, thermo, transData);

    //- Density [kg/m^3]
    const scalar rho = thermo.rho(species, T);

    //- ZRot
    const scalar ZRot298 = transData.ZRot298(species);

    //- Calculate Zrot
    const scalar ZRot = ZRot298 * F(scalar(298), LJP) / F(T, LJP);

    //- Calculate A
    const scalar A = 5./2. - rho * Dkk / mu;

    //- Calculate B
    const scalar B = ZRot + 2/M_PI * (5./3. * CvRot / R + rho * Dkk / mu);

    //- Calculate fTrans
    const scalar fTrans = 5./2. * (1 - (2*CvRot * A) / (M_PI * CvTrans * B));

    //- Calculate fVib
    const scalar fVib = rho * Dkk / mu;

    //- Calculate fRot
    const scalar fRot = fVib * ( 1 + (2 * A) / (M_PI * B));

    //- Calculate thermal conductivity 
    const scalar thermalConductivity =
        mu/MW * (fTrans * CvTrans + fRot * CvRot + fVib * CvVib);
    
    return 0;
}*/


AFC::scalar AFC::TransportCalc::thermalConductivityPolynomial
(
    const scalar T,
    const scalarField& polyCoeffs
) const
{
    return
    (
        exp
        (
            polyCoeffs[3]
          + polyCoeffs[2] * log(T)
          + polyCoeffs[1] * pow(log(T), 2)
          + polyCoeffs[0] * pow(log(T), 3)
        )
    );
}


void AFC::TransportCalc::fitThermalConductivity
(
    TransportData& transData,
    const Thermo& thermo
) const
{
    //- Species from chemistry
    const wordList& species = transData.chemistrySpecies();

    //- Fit for all species
    forAll(species, s)
    {
        //- Temperature field 
        scalarField T;

        //- Lambda field
        scalarField lambda;

        //- n values for fitting from 300K - 4300K
        const unsigned int n{10};

        const scalar dT = scalar(4000)/n;

        //- Fill the fields with necessary data
        for (size_t i = 0; i < n; ++i)
        {
            //- Temperature value
            const scalar Tv = scalar(300) + i * dT;

            //- Add temperature
            T.push_back(Tv);

            //- Add viscosity
            lambda.push_back
            (
                log(thermalConductivity(s, Tv, thermo, transData))
            );
        }

        //- Fitting 
        //  This could be done more beatuiful with objects
        //
        //  Polynom: ln(lambda) = D ln(T)^3 + C ln(T)^2  + B ln(T) + A 
        //  f1(T) = ln(T)^3
        //  f2(T) = ln(T)^2
        //  f3(T) = ln(T)
        //  f4(T) = 1

        //- Now we need to build a matrix A that looks like 
        //          1             2             3             4          5
        //1 | [f1(T)f1(T)]  [f1(T)f2(T)]  [f1(T)f3(T)]  [f1(T)f4(T)] [yf1(T)] |
        //2 | [f2(T)f1(T)]  [f2(T)f2(T)]  [f2(T)f3(T)]  [f2(T)f4(T)] [yf1(T)] |
        //3 | [f3(T)f1(T)]  [f3(T)f2(T)]  [f3(T)f3(T)]  [f3(T)f4(T)] [yf1(T)] |
        //4 | [f4(T)f1(T)]  [f4(T)f2(T)]  [f4(T)f3(T)]  [f4(T)f4(T)] [yf1(T)] |
        //
        //  [] means gauss summation; y = lambda value
        //  Matrix is symmetric in case of the first 4 rows and cols
        scalar gS_11{0};
        scalar gS_12{0};
        scalar gS_13{0};
        scalar gS_14{0};
        scalar gS_15{0};

        scalar gS_22{0};
        scalar gS_23{0};
        scalar gS_24{0};
        scalar gS_25{0};

        scalar gS_33{0};
        scalar gS_34{0};
        scalar gS_35{0};

        scalar gS_44{0};
        scalar gS_45{0};

        //- a) Calculate the gauss summations
        for (size_t j = 0; j < n; ++j)
        {
            // Temperature value
            const scalar& Tv = T[j];
            const scalar& lambdav = lambda[j];

            gS_11 += pow(log(Tv), 3) * pow(log(Tv), 3);
            gS_12 += pow(log(Tv), 3) * pow(log(Tv), 2);
            gS_13 += pow(log(Tv), 3) * pow(log(Tv), 1);
            gS_14 += pow(log(Tv), 3) * 1;
            gS_15 += lambdav * pow(log(Tv), 3);

            gS_22 += pow(log(Tv), 2) * pow(log(Tv), 2);
            gS_23 += pow(log(Tv), 2) * pow(log(Tv), 1);
            gS_24 += pow(log(Tv), 2) * 1;
            gS_25 += lambdav * pow(log(Tv), 2);

            gS_33 += pow(log(Tv), 1) * pow(log(Tv), 1);
            gS_34 += pow(log(Tv), 1) * 1;
            gS_35 += lambdav * pow(log(Tv), 1);

            gS_44 += 1 * 1;
            gS_45 += lambdav * 1;
        }

        //- b) Assign the values to the matrix A
        Matrix A(4,4);

        A(0,0,gS_11);
        A(0,1,gS_12);
        A(0,2,gS_13);
        A(0,3,gS_14);

        //- A_21 = A_12
        A(1,0,gS_12);
        A(1,1,gS_22);
        A(1,2,gS_23);
        A(1,3,gS_24);

        //  A_31 = A_13
        //  A_32 = A_23
        A(2,0,gS_13);
        A(2,1,gS_23);
        A(2,2,gS_33);
        A(2,3,gS_34);

        //  A_41 = A_14
        //  A_42 = A_24
        //  A_43 = A_34
        A(3,0,gS_14);
        A(3,1,gS_24);
        A(3,2,gS_34);
        A(3,3,gS_44);

        //- Calculate the inverse
        Matrix AI = A.inverse();

        //- Init the solution vector
        Vector b(4);

        b(0, gS_15);
        b(1, gS_25);
        b(2, gS_35);
        b(3, gS_45);

        //- Calculate the polynomial coefficients
        //  x = b * A^-1
        const Vector x = b * AI;

        //- Save the values
        transData.thermalConductivityPolyCoeffs(s, x);
    }
}
    

// * * * * * * * Calculation Functions For Binary Diffusivity * * * * * * * *//

AFC::scalar AFC::TransportCalc::binaryDiffusivity
(
    const word species1,
    const word species2,
    const scalar T,
    const Thermo& thermo,
    const TransportData& transData,
    const word method
) const
{
    if (method == "ChapmanAndEnskog")
    {
        return binaryDiffusivityChapmanAndEnskog
        (
            species1,
            species2,
            T,
            thermo,
            transData
        );
    }
    else if (method == "Polynomial")
    {
        //- Polynomial coefficients
        const scalarField& polyCoeffs =
            transData.binaryDiffusivityPolyCoeffs(species1, species2);

        return binaryDiffusivityPolynomial(T, polyCoeffs);
    }
    else
    {
        ErrorMsg
        (
            "   Method " + method + " is not available for calculation of\n"
            "   binary diffusifity. You have to implement it yourself\n"
            "   or ask for help",
            __FILE__,
            __LINE__
        );

        //- For compiler
        return -1;
    }
}


AFC::scalar AFC::TransportCalc::binaryDiffusivityChapmanAndEnskog
(
    const word species1,
    const word species2,
    const scalar T,
    const Thermo& thermo,
    const TransportData& transData
) const
{
    //- Pressure [bar]
    const scalar& p = thermo.p()/1e5;

    //- Molecular weight [g/mol]
    const scalar& MW1 = thermo.MW(species1);
    const scalar& MW2 = thermo.MW(species2);

    //- Lennard-Jones potential well depth eps/kb [K]
    const scalar& LJP1 = transData.LJP(species1);
    const scalar& LJP2 = transData.LJP(species2);

    //- Lennard-Jones collision diameter [Angstroms]
    const scalar& sigma1 = transData.LJCD(species1);
    const scalar& sigma2 = transData.LJCD(species2);

    //- Combined parameters 
    const scalar MW12 = 2 * pow((1/MW1 + 1/MW2), -1);
    const scalar LJP12 = sqrt(LJP1 * LJP2);
    const scalar sigma12 = (sigma1 + sigma2) / 2;

    //- Reduced temperature Ts12
    const scalar Ts12 = T * pow(LJP12, -1);

    //- Calculate collision integral omega 22
    const scalar omegaColl = reducedCollisionIntegralOmega11(Ts12);

    //- Calculate binary diffusivity [cm^2/s]
    const scalar Dij =
        0.002662 * pow(T, scalar(1.5))
      / (p * sqrt(MW12) * pow(sigma12, 2) * omegaColl);

    //- Return the binary diffusivity Dij [m^2/s] 
    //  Factor 1m^2 = 100cm * 100cm = 10000
    //  D is in cm^2/s so it means that it diffuse in 1 s for a special
    //  area - hence we have to divide by 10000 to get [m^2/s]
    return Dij / scalar(10000);
}


/*AFC::scalar AFC::TransportCalc::binaryDiffusivityWarnatz
(
    const word species1,
    const word species2,
    const scalar T,
    const Thermo& thermo,
    const TransportData& transData
) const
{
    //- Pressure [bar]
    const scalar& p = thermo.p()/1e5;

    //- Molecular weight [g/mol]
    const scalar& MW1 = thermo.MW(species1);
    const scalar& MW2 = thermo.MW(species2);

    //- Reduced molecular weight [g/mol]
    const scalar RMW = (MW1 + MW2) / (2 * MW1 * MW2);

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
    //const scalar Dij =
    //     2.662e-5 * sqrt(pow(T, 3) * 333; 

}*/


/*AFC::scalar AFC::TransportCalc::binaryDiffusivityCC
(
    const word species1,
    const word species2,
    const scalar T,
    const Thermo& thermo,
    const TransportData& transData
) const
{
    //- Stefan Boltzmann constant [J/m^2/K^4]
    const scalar& kB = AFC::Constants::kB;
    
    //- Pressure
    const scalar& p = thermo.p();

    //- Molecular weight [g/mol]
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

    //- Calculate collision integral omega 11
    const scalar omegaColl = reducedCollisionIntegralOmega11(Ts12);

    return 0;
}*/


AFC::scalar AFC::TransportCalc::binaryDiffusivityPolynomial
(
    const scalar T,
    const scalarField& pC
) const
{
    return
    (
        exp
        (
            pC[3]
          + pC[2] * log(T)
          + pC[1] * pow(log(T), 2)
          + pC[0] * pow(log(T), 3)
        )
    );
}


void AFC::TransportCalc::fitBinaryDiffusivity
(
    TransportData& transData,
    const Thermo& thermo
) const
{
    //- Species from chemistry
    const wordList& species = transData.chemistrySpecies();

    //- Fit for all species (first binary species)
    forAll(species, s1)
    {
        //- Second binary species 
        //  TODO (do not use all species // multpile components) ij = ji
        forAll(species, s2)
        {
            //- Temperature field 
            scalarField T;

            //- Binary diffusivity field Dij
            scalarField Dij;

            //- n values for fitting from 300K - 4300K
            const unsigned int n{10};

            const scalar dT = scalar(4000)/n;

            //- Pressure
            const scalar& p = thermo.p();

            //- Fill the fields with necessary data
            for (size_t i = 0; i < n; ++i)
            {
                //- Temperature value
                const scalar Tv = scalar(300) + i * dT;

                //- Add temperature
                T.push_back(Tv);

                //- Add binary ln(diffusivity)
                Dij.push_back
                (
                    log(binaryDiffusivity(s1, s2, Tv, thermo, transData))
                );
            }

            //- Fitting 
            //  This could be done more beatuiful with objects
            //
            //  Polynom: ln(Dij*p) = D ln(T)^3 + C ln(T)^2  + B ln(T) + A 
            //  f1(T) = ln(T)^3
            //  f2(T) = ln(T)^2
            //  f3(T) = ln(T)
            //  f4(T) = 1

            //- Now we need to build a matrix A that looks like 
            //          1             2             3             4          5
            //1 | [f1(T)f1(T)]  [f1(T)f2(T)]  [f1(T)f3(T)]  [f1(T)f4(T)] [yf1(T)] |
            //2 | [f2(T)f1(T)]  [f2(T)f2(T)]  [f2(T)f3(T)]  [f2(T)f4(T)] [yf1(T)] |
            //3 | [f3(T)f1(T)]  [f3(T)f2(T)]  [f3(T)f3(T)]  [f3(T)f4(T)] [yf1(T)] |
            //4 | [f4(T)f1(T)]  [f4(T)f2(T)]  [f4(T)f3(T)]  [f4(T)f4(T)] [yf1(T)] |
            //
            //  [] means gauss summation; y = ln(Dij*p) value
            //  Matrix is symmetric in case of the first 4 rows and cols
            scalar gS_11{0};
            scalar gS_12{0};
            scalar gS_13{0};
            scalar gS_14{0};
            scalar gS_15{0};

            scalar gS_22{0};
            scalar gS_23{0};
            scalar gS_24{0};
            scalar gS_25{0};

            scalar gS_33{0};
            scalar gS_34{0};
            scalar gS_35{0};

            scalar gS_44{0};
            scalar gS_45{0};

            //- a) Calculate the gauss summations
            for (size_t j = 0; j < n; ++j)
            {
                // Temperature value
                const scalar& Tv = T[j];
                const scalar& Dijv = Dij[j];

                gS_11 += pow(log(Tv), 3) * pow(log(Tv), 3);
                gS_12 += pow(log(Tv), 3) * pow(log(Tv), 2);
                gS_13 += pow(log(Tv), 3) * pow(log(Tv), 1);
                gS_14 += pow(log(Tv), 3) * 1;
                gS_15 += Dijv * pow(log(Tv), 3);

                gS_22 += pow(log(Tv), 2) * pow(log(Tv), 2);
                gS_23 += pow(log(Tv), 2) * pow(log(Tv), 1);
                gS_24 += pow(log(Tv), 2) * 1;
                gS_25 += Dijv * pow(log(Tv), 2);

                gS_33 += pow(log(Tv), 1) * pow(log(Tv), 1);
                gS_34 += pow(log(Tv), 1) * 1;
                gS_35 += Dijv * pow(log(Tv), 1);

                gS_44 += 1 * 1;
                gS_45 += Dijv * 1;
            }

            //- b) Assign the values to the matrix A
            Matrix A(4,4);

            A(0,0,gS_11);
            A(0,1,gS_12);
            A(0,2,gS_13);
            A(0,3,gS_14);

            //- A_21 = A_12
            A(1,0,gS_12);
            A(1,1,gS_22);
            A(1,2,gS_23);
            A(1,3,gS_24);

            //  A_31 = A_13
            //  A_32 = A_23
            A(2,0,gS_13);
            A(2,1,gS_23);
            A(2,2,gS_33);
            A(2,3,gS_34);

            //  A_41 = A_14
            //  A_42 = A_24
            //  A_43 = A_34
            A(3,0,gS_14);
            A(3,1,gS_24);
            A(3,2,gS_34);
            A(3,3,gS_44);

            //- Calculate the inverse
            Matrix AI = A.inverse();

            //- Init the solution vector
            Vector b(4);

            b(0, gS_15);
            b(1, gS_25);
            b(2, gS_35);
            b(3, gS_45);

            //- Calculate the polynomial coefficients
            //  x = b * A^-1
            const Vector x = b * AI;

            //- Save the values
            transData.binaryDiffusivityPolyCoeffs(s1, s2, x);
        }
    }
}
    

// * * * * * * * * * * * * Additional Functions * * * * * * * * * * * * * * *//

AFC::scalar AFC::TransportCalc::F(const scalar T, const scalar LJP) const
{
    //- LJP Lennard-Jones potential well depth epsilon/kB [K]
    //- Return the temperature modifier
    return
    (
        1
      + pow(M_PI, 1.5) / 2 * pow(LJP / T, 0.5)
      + (pow(M_PI, 2) / 4 + 2) * (LJP / T )
      + pow(M_PI, 1.5) * pow(LJP / T, 1.5)
    );
}

// ************************************************************************* //
