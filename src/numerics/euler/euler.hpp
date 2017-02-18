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
    AFC::Euler
    
Description
    Abstract AFC::Euler class for building and calculating matrices

SourceFiles
    numerics.cpp

\*---------------------------------------------------------------------------*/

#ifndef Euler_hpp
#define Euler_hpp

#include "typedef.hpp"
#include "matrix.hpp"
#include "jacobian.hpp"
#include "stepStatus.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class Euler Declaration
\*---------------------------------------------------------------------------*/

class Euler
{
    private:

        // Debug switch
        bool debug_{false};

        //- Pointer to Jacobian object
        Jacobian* jac_;

        //- Rate of concentration change
        mutable map<word, scalar> dcdt_;

        //- Jacobian Matrix
        mutable Matrix dcdc_;

        /*size_t n_{10};

        //- Euler constants and variables
        const size_t kMax_ = 12;
        const size_t iMax_ = 13;
        
        size_t kTarget_;

        List<size_t> pivotIndices_;
        List<size_t> nSeq_;

        scalar stepFactor_1;
        scalar stepFactor_2;
        scalar stepFactor_3;
        scalar stepFactor_4;
        scalar stepFactor_5;
        scalar kFactor1_;
        scalar kFactor2_;
        mutable scalar jacRedo_;
        mutable scalar theta_;

        mutable scalarField dcdt_;
        scalarField cpu_;
        mutable scalarField dtOpt_;
        mutable scalarField temp_;
        scalarField cSequence_;
        scalarField scale_;

        scalarField dc_;
        scalarField cTemp_;

        Matrix coeff_;
        mutable Matrix dcdc_;
        Matrix table_;

        mutable map<word, scalar> c0_;
        */


    public:

        //- Constructor 
        Euler
        (
            const Chemistry&
        );

        //- Destructor
        ~Euler();


        // Member functions

            //- Solve using euler algorithm
            void solve
            (
                const scalar,
                const scalar,
                map<word, scalar>&,
                scalar&,
                StepStatus&
            ) const;


    private:

        // Private member functions

            //- Comutes the j-th line of the extrapolation table
            bool seul
            (
                const scalar,
                const scalarField&,
                const scalar,
                const size_t,
                scalarField&,
                const scalarField&
            ) const;

            //- Polynomial extrapolation
            void extrapolate
            (
                const size_t,
                Matrix&,
                scalarField&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // Euler_hpp included

// ************************************************************************* //
