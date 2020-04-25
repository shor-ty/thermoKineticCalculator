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
    TKC::StepStatus
    
Description
    Abstract TKC::StepStatus class for building and calculating matrices

SourceFiles
    numerics.cpp

\*---------------------------------------------------------------------------*/

#ifndef StepStatus_hpp
#define StepStatus_hpp

#include "definitions.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TKC
{

/*---------------------------------------------------------------------------*\
                            Class StepStatus Declaration
\*---------------------------------------------------------------------------*/

class StepStatus
{
    public:
            
            //- Bool if we have to calculate or not
            bool forward_{false};

            //- Time step that we try
            scalar dtTry_{0};

            //- Time step that we did
            scalar dtDid_{0};

            //- First iteration?
            bool firstIter_{true};

            //- Last iteration?
            bool lastIter_{false};

            //- Reject 
            bool reject_{false};

            //- Previous reject 
            bool prevReject_{false};


    public:
        
        //- Constructor
        StepStatus(const scalar);

        //- Destructor
        ~StepStatus();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TKC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // StepStatus_hpp included

// ************************************************************************* //
