/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
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

#include "transport.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Transport::Transport(const string fileName, const Thermo& thermo)
:
    TransportCalc(fileName, thermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Transport::~Transport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AFC::Transport::fitCurves()
{
    Info<< " c-o Fit the transport properties to the polynomials\n" << endl;

    //- Fit the viscosity
    fitViscosity();

    //- Calculate thermal conductivity using gas kinetics
    fitThermalConductivity();
    
    //- Calculate binary diffusivity using gas kinetics
    fitBinaryDiffusivity();
}


// * * * * * * * * * * * * * * * Summary Functions * * * * * * * * * * * * * //

void AFC::Transport::summary(ostream& data) const
{
    data<< Header() << "\n";

    data<< " c-o Transport summary:\n"
        << "=======================\n\n"
        << "======================================================"
        << "===================================================\n"
        << "       Species      |   Geo. Config.    eps/kb    "
        << "    sigma          mu        alpha      zRot293    |\n"
        << "======================================================"
        << "===================================================\n";

    //  TODO switch to plot all or onyl the needed
    //- Just use the species that are used in chemistry
    const wordList& species_ = chemistrySpecies();

    data<<std::fixed;
    data.precision(4);

    forAll(species_, s)
    {
        data<< "  " << std::left << std::setw(15) << s << "   |"
            << std::right
            << "  " << std::setw(8) << geometricalConfig(s)
            << "  " << std::setw(13) << LJCD(s)
            << "  " << std::setw(12) << LJP(s)
            << "  " << std::setw(10) << muk(s)
            << "  " << std::setw(10) << alpha(s)
            << "  " << std::setw(10) << ZRot298(s)
            << "     |\n";
    }

    data<< "======================================================"
        << "===================================================\n";

    data<< "\n\n";

    forAll(species_, s)
    {
        data<< "============================================================="
            << "==\n"
            << " c-o Transport analyses for " << s << "\n"
            << "============================================================="
            << "==\n";

        data<< "\n\n-----------------------------------------------------";

        forAll(species_, ss)
        {
            if (s != ss)
            {
                data<< "-----------------";
            }
        }

        data<< "\n"
            << "      T    |         mu           lambda ";

        //- Species binary diffusion ss := second species
        forAll(species_, ss)
        {
            if (s != ss)
            {
                data<< "         Dij     ";
            }
        }

        data<< "\n           |" << std::setw(26) << " ";

        forAll(species_, ss)
        {
            if (s != ss)
            {
                const word Dij = s + "-" + ss;

                data<< std::setw(17) << Dij;
            }
        }

        data<< "\n";


        data<< "     [K]   |      [kg/m/s]       [W/m/K]    ";

        forAll(species_, ss)
        {
            if (s != ss)
            {
                data<< "     [m^2/s]     ";
            }
        }

        data<< "\n";

        data<< "---------------------------------------------------"
            << "--";

        forAll(species_, ss)
        {
            if (s != ss)
            {
                data<< "-----------------";
            }
        }

        data<< "\n";

        data<<std::scientific;

        for(int i=300; i<=3000; i+=100)
        {
            data<< "  " << std::setw(6) << i << "   |" 
                << "  " << std::setw(13) << viscosity(s, i) 
                << "  " << std::setw(13) << thermalConductivity(s, i);

            //- Species binary diffusion ss := second species
            forAll(species_, ss)
            {
                if (s != ss)
                {
                    data<< "  " << std::setw(13) << binaryDiffusivity(s, ss, i)
                        << "  ";
                }
            }

            data<< "\n";
        }

        data<< "---------------------------------------------------";

        forAll(species_, ss)
        {
            if (s != ss)
            {
                data<< "-----------------";
            }
        }

        data<< "\n\n\n\n";
    }
}


void AFC::Transport::summaryFittingProcedure(ostream& data) const
{
    data<< Header() << "\n";

    data<< " c-o Fitting procedure summary:\n"
        << "===============================\n\n"
        << "======================================================\n"
        << "   The fitting polynomials that are used\n"
        << "======================================================\n\n"
        << "   Viscosity            : "
        << " mu = e^( D ln(T)^3 + C ln(T)^2 + B ln(T) + A )      [kg/m/s]\n\n"
        << "   Thermal Conductivity : "
        << " lambda = e^( D ln(T)^3 + C ln(T)^2 + B ln(T) + A )  [W/m/K]\n\n"
        << "   Binary Diffusivity   : "
        << " Dij = e^( D ln(T)^3 + C ln(T)^2 + B ln(T) + A ) / p [m^2/s]\n\n"
        << "\n\n======================================================\n"
        << "   Viscosity polynomial coefficients\n"
        << "======================================================\n\n";

    //- Table of viscosity polynomial coefficients
    data<< "------------------------------------------------------------"
        << "------------------------------------------------------------\n"
        << "  Species        |        A               B               C"
        << "               D             mu(298K)        mu(1000K)   |\n"
        << "------------------------------------------------------------"
        << "------------------------------------------------------------\n";
    
    const wordList& species = chemistrySpecies();

    forAll(species, s)
    {
        //- Polynomial Coefficients
        const scalarField& pC = viscosityPolyCoeffs(s);

        data<< "  " << std::left << std::setw(15)<< s << "|" << std::right
            << std::setw(16) << pC[3]
            << std::setw(16) << pC[2]
            << std::setw(16) << pC[1]
            << std::setw(16) << pC[0]
            << std::setw(16) << viscosity(s, 298, "Polynomial")
            << std::setw(16) << viscosity(s, 1000, "Polynomial") 
            << "  |\n";
    }

    data<< "------------------------------------------------------------"
        << "------------------------------------------------------------\n\n"
        << "\n\n======================================================\n"
        << "   Thermal conductivity coefficients\n"
        << "======================================================\n\n";

    //- Table of thermal conductivity polynomial coefficients
    data<< "------------------------------------------------------------"
        << "------------------------------------------------------------\n"
        << "  Species        |        A               B               C"
        << "               D           lambda(298K)    lambda(1000K) |\n"
        << "------------------------------------------------------------"
        << "------------------------------------------------------------\n";
    
    forAll(species, s)
    {
        //- Polynomial Coefficients
        const scalarField& pC = thermalConductivityPolyCoeffs(s);

        data<< "  " << std::left << std::setw(15)<< s << "|" << std::right
            << std::setw(16) << pC[3]
            << std::setw(16) << pC[2]
            << std::setw(16) << pC[1]
            << std::setw(16) << pC[0]
            << std::setw(16) << thermalConductivity(s, 298, "Polynomial")
            << std::setw(16) << thermalConductivity(s, 1000, "Polynomial") 
            << "  |\n";
    }

    data<< "------------------------------------------------------------"
        << "------------------------------------------------------------\n\n"
        << "\n\n======================================================\n"
        << "   Binary diffusion coefficients\n"
        << "======================================================\n\n";

    //- Table of binary diffusion polynomial coefficients
    data<< "------------------------------------------------------------"
        << "------------------------------------------------------------"
        << "---------------------\n"
        << "  Species 1         Species 2         |        A               B"
        << "               C"
        << "               D              Dij(298K)       Dij(1000K) |\n"
        << "------------------------------------------------------------"
        << "------------------------------------------------------------"
        << "---------------------\n";
    
    //- First species
    forAll(species, s1)
    {
        //- Second species
        forAll(species, s2)
        {
            //- Polynomial Coefficients
            const scalarField& pC = binaryDiffusivityPolyCoeffs(s1, s2);

            data<< "  " << std::left << std::setw(18)<< s1
                << std::setw(18) << s2 << "|" << std::right
                << std::setw(16) << pC[3]
                << std::setw(16) << pC[2]
                << std::setw(16) << pC[1]
                << std::setw(16) << pC[0]
                << std::setw(16)
                << binaryDiffusivity(s1, s2, 298, "Polynomial")
                << std::setw(16)
                << binaryDiffusivity(s1, s2, 1000, "Polynomial") 
                << "  |\n";
        }
        data<< "------------------------------------------------------------"
            << "------------------------------------------------------------"
            << "---------------------\n";
    }

    data<< "-------------------------------------------------------------"
        << "-------------------------------------------------------------\n";


}


// ************************************************************************* //