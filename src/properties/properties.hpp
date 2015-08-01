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
    AFC::Properties
    
Description
    Abstract AFC::Properties class for properties data and calculation

SourceFiles
    properties.cpp

\*---------------------------------------------------------------------------*/

#ifndef Properties_hpp
#define Properties_hpp

#include "propertiesReader.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Properties Declaration
\*---------------------------------------------------------------------------*/

class Properties
{
    private:

        // Private Data
        
            //- Mixture fraction discrete points
            int mfPoints_{0};

            //- Varianz of mixture fraction discrete points
            int vmfPoints_{0};

            //- ScalarDissipation rate
            scalarField sDR_;

            //- Oxidizer species
            wordList speciesOxidizer_;

            //- Composition of oxidizer [mol]
            map<word, scalar> oxidizerCompMol_;

            //- Fuel species
            wordList speciesFuel_;

            //- Composition of fuel [mol]
            map<word, scalar> fuelCompMol_;

            //- Temperature of oxidizer stream
            scalar TOxidizer_{0};

            //- Temperature of fuel stream
            scalar TFuel_{0};


        // Debug
        const bool debug{false};


    public:

        //- Constructor
        Properties
        (
            const string& 
        );

        //- Destructor
        ~Properties();


        // Member functions
        
            //- Insert mfPoints_
            void insertMFPoints
            (
                const int&
            );
        
            //- Insert vmfPoints_
            void insertVMFPoints
            (
                const int&
            );

            //- Insert scalar dissipation rates
            void insertScalarDissipationRates
            (
                const scalar&
            );

            //- Insert temperature oxidizer
            void insertTemperatureOxidizer
            (
                const scalar&
            );

            //- Insert temperature fuel
            void insertTemperatureFuel
            (
                const scalar&
            );
            
            //- Insert oxidizer mol composition 
            void insertCompositionOxidizerMol
            (
                const word&,
                const scalar&
            );
            
            //- Insert fuel mol composition 
            void insertCompositionFuelMol
            (
                const word&,
                const scalar&
            );
            

        // Check functions

            //- Check if all data are set
            void check();


        // Return functions

            //- Return oxidizer species
            wordList speciesOxidizer() const;

            //- Return fuel species
            wordList speciesFuel() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Properties_hpp included

// ************************************************************************* //
