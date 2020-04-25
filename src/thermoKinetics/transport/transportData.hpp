/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     T hermo       | Open Source Thermo-Kinetic Library
  c-o-o-c     K iknetic     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Constructor.

    TKC is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or 
    (at your option) any later version.

    TKC is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TKC; if not, see <http://www.gnu.org/licenses/>

Class
   TKC::TransportData
   
Description
    This class contains all transport data 

SourceFiles
    transportData.cpp

\*---------------------------------------------------------------------------*/

#ifndef TransportData_hpp
#define TransportData_hpp

#include "definitions.hpp"
#include "vector.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*--------------------------------------------------------------------------*\
                         Class TransportData Declaration
\*--------------------------------------------------------------------------*/

class TransportData
{
    private:

        // Private data

            //- List of species of transport file
            wordList species_;

            //- List of species used in the elementar reactions
            wordList chemicalFormula_;

            //- List of species used for elementar reactions
            wordList chemistrySpecies_;

            //- List for geometrical configuration
            //  Definiton:
            //      + 0 molecule is single atom
            //      + 1 molecule is linear
            //      + 2 molecule is nonlinear
            map<word, int> geoConfig_;

            //- List for Lennard-Jones potential well depth eps/kb [K]
            map<word, scalar> lenJonPot_;

            //- List for Lennard-Jones collision dimater C in Angstroms
            //  Also denoted by sigma 
            map<word, scalar> lenJonCollDia_;

            //- List for dipole moment mu in debey
            //  a debey is 10^-18 [cm^3/2]
            map<word, scalar> dipMom_;

            //- List for polarizability alpha in cubic Angstroms
            map<word, scalar> alpha_;

            //- List for rotational relaxation collision number Zrot
            //  at 298K
            map<word, scalar> ZRot298_;

            //- List of species binary combinations
            map<word, List<word> > binarySpeciesCombinations_;


        // Transport data for fitting procedure 

            //- Pure species viscosity
            //  Species, Coeffs A-D, Value
            map<word, scalarField> viscosity_;

            //- Pure species thermal conductivity 
            //  Species, Coeffs A-D, Value
            map<word, scalarField> thermalConductivity_;

            //- Pure binary diffusion coefficients
            //  Species1 <- ID from species_, Species2, Coeffs A-D, Value
            List<map<word, scalarField> > binaryDiffusivity_;


    public:

        //- Constructor 
        TransportData(const string);

        //- Destructor
        ~TransportData();


        // Member functions

            //- Fitting polynomials coefficients for the viscosity
            void viscosityPolyCoeffs(const word, const Vector&);

            //- Return the polynomial coefficients for the viscosity
            scalarField viscosityPolyCoeffs(const word) const;

            //- Fitting polynomials coefficients for the thermal conductivity 
            void thermalConductivityPolyCoeffs(const word, const Vector&);

            //- Return the polynomial coefficients for the thermal conductivity 
            scalarField thermalConductivityPolyCoeffs(const word) const;

            //- Fitting polynomials coefficients for the binary diffusivity 
            void binaryDiffusivityPolyCoeffs
            (
                const word,
                const word,
                const Vector&
            );

            //- Return the polynomial coefficients for the binary diffusivity 
            scalarField binaryDiffusivityPolyCoeffs
            (
                const word,
                const word
            ) const;


        // Insert functions, from TransportReader:: delegated 

            //- Insert species
            void insertSpecies(const word);

            //- Insert species used in elementar reactions
            void insertChemistrySpecies(const wordList&);

            //- Insert geometrical configuration
            void insertGeoConfig(const int);

            //- Insert Lennard-Jones potential
            void insertLenJonPot(const scalar);

            //- Insert Lennard-Jones collision diameter
            void insertLenJonCollDia(const scalar);

            //- Insert dipole momentum
            void insertDipMom(const scalar);

            //- Insert polarizability alpha
            void insertAlpha(const scalar);

            //- Insert rotational relaxation collision number
            void insertZRot298(const scalar);

            //- Insert binary species combination
            void insertBinarySpeciesCombinations(const word, const word);


        // Insert function from afc.cpp

            //- Insert chemical species
            void chemicalFormula(const wordList&);


        // Update and Manipulation Functions

            //- Build the binary species combinations
            void binarySpeciesCombinations();


        // Return functions
        
            //- Return all species which are stored in the transport object
            wordList species() const;

            //- Return chemical formula of all stored species
            wordList chemicalFormula() const;

            //- Return chemical formula of species s
            word chemicalFormula(const word) const;

            //- Return all species which are used in elementar reactions
            wordList chemistrySpecies() const;


            //- Return the geometrical configuration
            int geometricalConfig(const word) const;

            //- Return Lennard-Jones collision diameter [Angstroms]
            scalar LJCD(const word) const;

            //- Return Lennard-Jones potential well depth eps/kb [K]
            scalar LJP(const word) const;

            //- Return the dipole moment [debey]
            scalar muk(const word) const;

            //- Return polarizability alpha [cubic Angstroms]
            scalar alpha(const word) const;

            //- Return rotational relaxation collision number Zrot
            scalar ZRot298(const word) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // TransportData_hpp included

// ************************************************************************* //
