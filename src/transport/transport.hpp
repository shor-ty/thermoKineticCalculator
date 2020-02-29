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

Class
    AFC::Transport
    
Description
    Abstract AFC::Transport class for transport data and calculation

SourceFiles
    transport.cpp

\*---------------------------------------------------------------------------*/

#ifndef Transport_hpp
#define Transport_hpp

#include "transportReader.hpp"
#include "transportData.hpp"
#include "transportCalc.hpp"
#include "thermo.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Transport Declaration
\*---------------------------------------------------------------------------*/

class Transport
{
    private:

        // Private reference data

            //- Transport Data
            TransportData transData_;

            //- Transport Calc
            TransportCalc transCalc_;
        
            //- Thermodynamic Object 
            const Thermo& thermo_;


    public:

        //- Constructor
        Transport(const string, const Thermo&);

        //- Destructor
        ~Transport();


        // Calculation Functions

            //- Insert chemistry species to transportData
            void insertChemistrySpecies(const wordList&);
        
            //- Calculate fitting coefficients for polynomials
            void fitCurves();

            //- Calculate viscosity of species s [Pa s]
            scalar viscosity 
            (
                const word,
                const scalar,
                const word method = "Hirschfelder"
            ) const;

            //- Calculate thermal conductivity of species s [W/m/K]
            scalar thermalConductivity 
            (
                const word,
                const scalar,
                const word method = "Warnatz"
            ) const;

            //- Calculate the binary diffusivity coefficient Dij of species ij
            scalar binaryDiffusivity
            (
                const word,
                const word,
                const scalar,
                const word method = "ChapmanAndEnskog"
            ) const;

        // Return functions

            //- Return species as wordList
            wordList species() const;

            //- Return the species used in elementar reactions
            wordList chemistrySpecies() const;

            //- Return chemical formula for species as wordList
            wordList chemicalFormula() const;

            //- Return chemical formula for species s
            word chemicalFormula(const word) const;

            //- Return geometrical configuration of species s
            int geometricalConfig(const word) const;

            //- Return Lennard-Jones collision diameter [Angstroms]
            scalar LJCD(const word) const;

            //- Return Lennard-Jones potential well depth eps/kb [K]
            scalar LJP(const word) const;

            //- Return the dipole moment [debey]
            scalar muk(const word) const;

            //- Return polarizability alpha [cubic Angstroms]
            scalar alpha(const word) const;

            //- Return rotational relaxation collision number Zrot at 298 K
            scalar ZRot298(const word) const;

            //- Return the viscosity polynomial coefficients
            scalarField viscosityPolyCoeffs(const word) const;

            //- Return the thermal conductivity polynomial coefficients
            scalarField thermalConductivityPolyCoeffs(const word) const;

            //- Return the binary diffusivity polynomial coefficients
            scalarField binaryDiffusivityPolyCoeffs
            (
                const word,
                const word
            ) const;


        // Summary functions

            //- Build the summary
            void summary(ostream&) const;

            //- Build the fitting procedure summary
            void summaryFittingProcedure(ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Transport_hpp included

// ************************************************************************* //
