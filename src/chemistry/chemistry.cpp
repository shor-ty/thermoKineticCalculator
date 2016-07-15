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

#include "chemistry.hpp"
#include "thermo.hpp"
#include "constants.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Chemistry::Chemistry
(
    const string& fileName 
)
{
    if (debug_)
    {
        Info<< "Constructor Chemistry\n" << endl;
    }

    ChemistryReader chemReader(fileName);
    
    chemReader.read(chemData_);

    //- Reactions that include species i
    createSpeciesInReaction();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Chemistry::~Chemistry()
{
    if (debug_)
    {
        Info<< "Destructor Chemistry\n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool AFC::Chemistry::thermo()
{
    return (chemData_.thermo());
}


// * * * * * * * * * * * * * Calculation Functions * * * * * * * * * * * * * //

AFC::scalar AFC::Chemistry::calculateOmega
(
    const word& species,
    const scalar& T,
    map<word, scalar>& con,
    const Thermo& thermo
)
{
    //- Calculate source term omega
    //return chemCalc_.calculateOmega(species, T, con, thermo, chemData_);
    return 0;
}


void AFC::Chemistry::calculateKf
(
    const int& r,
    const scalar& T
)
{
    //chemCalc_.calculateKf(r, T, chemData_);
}


void AFC::Chemistry::calculateKc
(
    const int& r,
    const scalar& T,
    const Thermo& thermo
)
{
    //chemCalc_.calculateKc(r, T, thermo, chemData_);
}


void AFC::Chemistry::calculateKb()
{
    //chemCalc_.calculateKb(chemData_);
}


// * * * * * * * * * * * * * * * Update Functions  * * * * * * * * * * * * * //

/*void AFC::Chemistry::updateM
(
    const scalar& M
)
{
    chemData_.updateM(M);
}*/


// * * * * * * * * * * * * * * * Create Functions  * * * * * * * * * * * * * //

void AFC::Chemistry::createSpeciesInReaction()
{
    if (debug_)
    {
        Info<< " --> AFC::Chemistry::createSpeciesInReaction" << endl;
    }

    //- Species list
    const wordList& species = chemData_.species();

    //- Reaction no.
    const int& nReac = chemData_.nReac();

    //- List<wordList> that contains all species in each reaction
    const List<wordList>& speciesInReaction = chemData_.speciesInReaction();

    forAll(species, s)
    {
        //- If found species in reaction -> true
        bool found{false};

        //- Loop through all elementar reactions 
        for(int r=0; r<nReac; r++)
        {
            //- Loop through the species i in elementar reaction r
            for(unsigned int i=0; i<speciesInReaction[r].size(); i++)
            {
                if
                (
                  ! found
                 && speciesInReaction[r][i] == s
                )
                {
                    found = true;

                    //- Insert reaction no to matrix
                    chemData_.setReacNumbers(s, r);
                }
            }

            //- Reset
            found = false;
        }
    }
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

bool AFC::Chemistry::BR
(
    const int& reacNo
) const
{
    return chemData_.BR(reacNo);
}


bool AFC::Chemistry::TBR
(
    const int& reacNo
) const
{
    return chemData_.TBR(reacNo);
}


bool AFC::Chemistry::LOW
(
    const int& reacNo
) const
{
    return chemData_.LOW(reacNo);
}


bool AFC::Chemistry::TROE
(
    const int& reacNo
) const
{
    return chemData_.TROE(reacNo);
}


bool AFC::Chemistry::SRI
(
    const int& reacNo
) const
{
    return chemData_.SRI(reacNo);
}


bool AFC::Chemistry::ENHANCED
(
    const int& reacNo
) const
{
    return chemData_.ENHANCED(reacNo);
}


AFC::wordList AFC::Chemistry::species() const
{
    return chemData_.species();
}


AFC::wordList AFC::Chemistry::elements() const
{
    return chemData_.elements();
}



/*

int AFC::Chemistry::nReac() const
{
    return chemData_.nReac();
}
*/


AFC::string AFC::Chemistry::elementarReaction
(
    const int& r
) const
{
    return chemData_.elementarReaction(r);
}


AFC::List<AFC::string> AFC::Chemistry::elementarReaction() const
{
    return chemData_.elementarReaction();
}


/*AFC::scalarField AFC::Chemistry::reacNoForSpecies
(
    const int& s
) const
{
    return chemData_.reacNoForSpecies(s);
}


AFC::scalarField AFC::Chemistry::k() const
{
    return chemData_.k();
}*/


/*AFC::wordMatrix AFC::Chemistry::speciesInReactions() const
{
    return chemData_.speciesInReactions();
}*/


/*AFC::scalar AFC::Chemistry::k
(
    const int& reacNo 
) const
{
    return chemData_.k(reacNo);
}*/


// * * * * * * * * * * * * * * * Summary Functions * * * * * * * * * * * * * //

void AFC::Chemistry::summary
(
    ostream& data 
) const
{

    const wordList& elements = chemData_.elements();
    const wordList& species = chemData_.species();
    const wordList& reactions = chemData_.elementarReaction();

    //- Get amount of LOW TROE SRI ENHANCED BR reactions
    unsigned int LOW{0};
    unsigned int TROE{0};
    unsigned int SRI{0};
    unsigned int ENH{0};
    unsigned int BR{0};
    unsigned int TBR{0};

    forEach(reactions, r)
    {
        if (chemData_.LOW(r))
        {
            LOW++;
        }
        if (chemData_.TROE(r))
        {
            TROE++;
        }
        if (chemData_.SRI(r))
        {
            SRI++;
        }
        if (chemData_.ENHANCED(r))
        {
            ENH++;
        }
        if (chemData_.BR(r))
        {
            BR++;
        }
        if (chemData_.TBR(r))
        {
            TBR++;
        }
    }

    //- General stuff
    data<< " c-o Reaction summary:\n";
    data<< " =====================\n\n";
    data<< " Number of elements             " << elements.size()  << "\n";
    data<< " Number of species              " << species.size()   << "\n";
    data<< " Number of elementar reactions  " << reactions.size() << "\n";
    data<< "   |\n";
    data<< "   |--> Number of backward reactions   " << BR   << "\n";
    data<< "   |--> Number of ThirdBody reactions  " << TBR  << "\n";
    data<< "   |--> Number of ENHANCED reactions   " << ENH  << "\n";
    data<< "   |--> Number of LOW reactions        " << LOW  << "\n";
    data<< "   |--> Number of TROE reactions       " << TROE << "\n";
    data<< "   |--> Number of SRI reactions        " << SRI  << "\n";

    data<< "\n\n Elements used: \n   | \n";

    forEach(elements, e)
    {
        data<< "   |-->  " << elements[e] << "\n";
    }

    data<< "\n\n Species used:\n   | \n";

    forEach(species, s)
    {
        data<< "   |-->  " << species[s] << "\n";
    }
    
    data<< "\n\n Reactions used: \n   | \n";

    forEach(reactions, r)
    {
        data<< "   |-->  (" << r+1 << ")  " << reactions[r] << "\n";
    }

    data<< "\n\n More detailed kinetics for reactions\n\n";

    forEach(reactions, r)
    {
        data<< "   Reaction " << r+1 << ":  " << reactions[r] << "\n";
        data<< "   =============================================\n\n";
        data<< "      A:  " << chemData_.arrheniusCoeffs(r)[0] << " ";
        data<< "\n";
        data<< "      n:  " << chemData_.arrheniusCoeffs(r)[1] << " ";
        data<< "\n";
        data<< "      Ea: " << chemData_.arrheniusCoeffs(r)[2] << " ";
        data<< "cal/mol\n\n";

        if(chemData_.TBR(r))
        {

            if(chemData_.LOW(r))
            {
                data<< "     LOW Coeffs (for high pressure)\n\n";

                const List<scalar>& LOWCoeffs = chemData_.LOWCoeffs(r);

                data<< "         A:  " << LOWCoeffs[0] << " ";
                data<< "\n";
                data<< "         n:  " << LOWCoeffs[1] << " ";
                data<< "\n";
                data<< "         Ea: " << LOWCoeffs[2] << " ";
                data<< "cal/mol\n\n";
            }

            if(chemData_.TROE(r))
            {
                data<< "     TROE Coeffs\n\n"; 

                const List<scalar>& TROECoeffs = chemData_.TROECoeffs(r);

                data<< "        a:  " << TROECoeffs[0] << "\n"; 
                data<< "        b:  " << TROECoeffs[1] << "\n"; 
                data<< "        c:  " << TROECoeffs[2] << "\n"; 
                data<< "        d:  " << TROECoeffs[3] << "\n"; 
            }
        }

        data<< "\n";
    }
}


// ************************************************************************* //
