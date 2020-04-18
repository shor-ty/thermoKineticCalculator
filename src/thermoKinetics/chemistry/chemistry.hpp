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
    AFC::Chemistry

Description
    Abstract AFC::Chemistry class for chemistry data and calculation

SourceFiles
    chemistry.cpp

\*---------------------------------------------------------------------------*/

#ifndef Chemistry_hpp
#define Chemistry_hpp

#include "chemistryCalc.hpp"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Chemistry Declaration
\*---------------------------------------------------------------------------*/

class Chemistry
:
    public ChemistryCalc
{

    public:

        //- Constructor
        Chemistry(const string, const Thermo&);

        //- Destructor
        ~Chemistry();


        // Member Functions

            //- Check if all species that are used in the chemistry are
            //  available in the thermo object
            void checkSpecies() const;


        // Calculation Functions


            //- Calculate the source term of all species (omega) and return it
            //  omega = dcdt
            map<word, scalar> omega
            (
                const scalar,
                const map<word, scalar>&
            ) const;



        // Update Functions


        // Create Functions

            //- Create field that contains all reaction no. for each species
            void buildSpeciesInReactionTable();

            //- Insert reaction no. for species delegated to CHEMISTRYDATA
            //void insertReacNoForSpecies(const int, const int);


        // Return Functions



        //- Summary Functions

            //- Build the summary as ostream
            void summary(ostream&) const;

            //- Build chemical table for reaction r (for summary)
            void chemicalTable(ostream&) const;

            //- Build table for kf and kb based on general coefficients
            //  or if we specify "LOW" then based on low coefficients
            void buildTablekf(const int, ostream&,const bool LOW = false) const;

            //- Build TROE table (calculate F_cent and logF)
            void buildTROETable(const int, ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Chemistry_hpp included

// ************************************************************************* //
