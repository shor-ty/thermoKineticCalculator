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

#include "transportReader.hpp" 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::TransportReader::TransportReader
(
    const string& file
)
:
    file_(file)

{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::TransportReader::~TransportReader()
{
    Info<< "Destructor TransportReader\n";
}


// * * * * * * * * * * * * * Runtime object creator  * * * * * * * * * * * * //

void AFC::TransportReader::newTransportData()
{
    pTrD_ = smartPtr<TransportData>(new TransportData());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

AFC::smartPtr<AFC::TransportData> AFC::TransportReader::readTransport()
{
    Info<< " c-o Reading transport data\n" << endl;

    newTransportData();

    const auto fileContent = readFile(file_);

    //- Reading transport properties
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
            pTrD_->insertSpecies(tmp[0]);

            pTrD_->insertGeoConfig(stoi(tmp[1]));

            pTrD_->insertLenJonPot(stod(tmp[2]));

            pTrD_->insertLenJonCollDia(stod(tmp[3]));

            pTrD_->insertDipMom(stod(tmp[4]));

            pTrD_->insertPol(stod(tmp[5]));

            pTrD_->insertRotRelCollNumb(stod(tmp[6]));
        }
    }

    return std::move(pTrD_);
}


// * * * * * * * * * * * * * * Helper functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * Data manipulation functions * * * * * * * * * * * //



// ************************************************************************* //
