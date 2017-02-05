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
    AFC::StepStatus
    
Description
    Abstract AFC::StepStatus class for building and calculating matrices

SourceFiles
    numerics.cpp

\*---------------------------------------------------------------------------*/

#ifndef StepStatus_hpp
#define StepStatus_hpp

#include "typedef.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace AFC
{

/*---------------------------------------------------------------------------*\
                            Class StepStatus Declaration
\*---------------------------------------------------------------------------*/

class StepStatus
{
    private:

        // Debug switch
        bool debug_{false};


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

} // End namespace AFC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // StepStatus_hpp included

// ************************************************************************* //
