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

\*---------------------------------------------------------------------------*/

#include "definitions.hpp"
#include "timeReader.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TKC::TimeReader::TimeReader(const string& file)
:
    file_(file)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TKC::TimeReader::~TimeReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void TKC::TimeReader::read(Time& time)
{
    Info<< " c-o Reading time control data\n"
        << "     >> " << file_ << "\n" << endl;

    const auto fileContent = readFile(file_);

    //- Reading file
    for (unsigned int line=0; line < fileContent.size(); line++)
    {
        string ttmp = fileContent[line];

        //- Remove any comment
        removeComment(ttmp);

        stringList tmp = splitStrAtWS(ttmp);

        //- If line is not empty, proceed
        if (!tmp.empty())
        {
            if (tmp[0] == "writeTime")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No value for writeTime is specified or it "
                        "is not a correct type (" + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                time.writeTime(stod(tmp[1]));
            }
            else if (tmp[0] == "initialKineticDeltaT")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No value for initialKineticDeltaT is specified or it "
                        "is not a correct type (" + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                time.dTKinetic(stod(tmp[1]));
                time.dTFlow(stod(tmp[1]));
            }
            else if (tmp[0] == "initialFlowDeltaT")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No value for initialFlowDeltaT is specified or it "
                        "is not a correct type (" + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                time.dTFlow(stod(tmp[1]));
            }
            else if (tmp[0] == "endTime")
            {
                if (tmp[1].empty())
                {
                    ErrorMsg
                    (
                        "No value for endTime is specified or it "
                        "is not a correct type (" + file_ + ")",
                        __FILE__,
                        __LINE__
                    );
                }

                time.endTime(stod(tmp[1]));
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
