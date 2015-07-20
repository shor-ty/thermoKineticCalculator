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

#include "mixtureFractionReader.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::MixtureFractionReader::MixtureFractionReader
(
    const string& file
)
:
    file_(file)

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::MixtureFractionReader::~MixtureFractionReader()
{
    Info<< "Destructor MixtureFractionReader\n" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AFC::MixtureFractionReader::read
(
    MixtureFraction& data
)
{
    Info<< " c-o Reading AFCDict\n" << endl;

    const auto fileContent = readFile(file_);

    //- Reading afcDict
    for (unsigned int line=0; line < fileContent.size(); line++)
    {
        stringList tmp = splitStrAtWS(fileContent[line]);

        //- If line is not empty and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            //- MixtureFraction points
            if (tmp[0] == "mixtureFractionPoints")
            {
                data.insertMFPoints(stoi(tmp[1]));
            }

        }
    }
}


// * * * * * * * * * * * * * * Helper functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * Data manipulation functions * * * * * * * * * * * //


// ************************************************************************* //
