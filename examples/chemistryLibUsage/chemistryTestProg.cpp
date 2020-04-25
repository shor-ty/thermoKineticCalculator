/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     T hermo       | Open Source Thermo-Kinetic Library
  c-o-o-c     K iknetic     |
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License
    This file is part of Automatic Flamelet Creator.

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

Description
    Example how to use the TKC modules for using the thermodynamic library
    Here we calculate the reaction enthalpy of different species out of the
    box as well as the free Gibbs energy


\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "definitions.hpp"
#include "chemistry.hpp"
#include "thermo.hpp"

using namespace TKC;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main ()
{


    const std::clock_t startTime = clock();

    Info<< Header() << endl;

    //- To use the chemistry we need a thermodynamic database first
    const Thermo thermo("files/NASA");

    //- Create the chemistry object
    const Chemistry chem("files/h2Kinetics", thermo);

    //-------------------------------------------------------------------------
    Info<< " c-o Number of elementar reactions: " << chem.nReac() << "\n"
        << " c-o Number of duplicated reactions: " << chem.nDuplicated() << "\n"
        << " c-o Number of ignored reactions: " << chem.nIgnored() << "\n"
        << endl;

    Info<< " c-o The chemistry contains the following elementar reactions:\n";
    const wordList& reactions = chem.elementarReaction();

    forAll(reactions, r)
    {
        Info<< " --> " << r << "\n";
    }
    Info<< endl;

    unsigned int nForwardReactions{0};
    unsigned int nBackwardReactions{0};
    unsigned int nBothReactionWays{0};
    forEach(reactions, i)
    {
        if (chem.FR(i) && chem.BR(i))
        {
            ++nBothReactionWays;
        }
        else if (chem.FR(i) && !chem.BR(i))
        {
            ++nForwardReactions;
        }
        else if (!chem.FR(i) && chem.BR(i))
        {
            ++nBackwardReactions;
        }
    }

    Info<< " c-o Number of irreversible forward reactions: "
        << nForwardReactions << "\n"
        << " c-o Number of irreversible backward reactions: "
        << nBackwardReactions << "\n"
        << " c-o Number of reversible reactions: " << nBothReactionWays << "\n"
        << endl;


    //-------------------------------------------------------------------------
    Info<< " c-o The species within reaction " << chem.elementarReaction(0)
        << " are:\n";

    //- Take care as 0 does mean integer and 0 pointer
    const wordList allSpeciesOfReaction = chem.speciesInReaction(0);

    forAll(allSpeciesOfReaction, s)
    {
        Info<< " --> " << s << "\n";
    }
    Info<< endl;

    //-------------------------------------------------------------------------
    Info<< " c-o The educt site of this reaction consists of:\n";
    const wordList& educts = chem.educts(0);

    forAll(educts, e)
    {
        Info<< " --> " << e << "\n";
    }
    Info<< endl;

    //-------------------------------------------------------------------------
    Info<< " c-o The product site of this reaction consists of:\n";
    const wordList& products = chem.products(0);

    forAll(products, p)
    {
        Info<< " --> " << p << "\n";
    }
    Info<< endl;

    //-------------------------------------------------------------------------
    Info<< " c-o The stochiometric factors are given as:\n";
    const map<word, int>& nuE = chem.nuEducts(int(0));
    const map<word, int>& nuP = chem.nuProducts(int(0));

    forAll(educts, e)
    {
        Info<< " --> " << e << " = " << nuE.at(e) << "\n";
    }
    forAll(products, p)
    {
        Info<< " --> " << p << " = " << nuP.at(p) << "\n";
    }
    Info<< endl;

    //- Reaction number
    unsigned int r = 0;

    //-------------------------------------------------------------------------
    Info<< " c-o This reaction is a ";
    if (chem.BR(r) && !chem.FR(r)) { Info<< "backward"; }
    if (chem.FR(r) && !chem.BR(r)) { Info<< "forward"; }
    if (chem.FR(r) && chem.BR(r)) { Info<< "reversible"; }
    Info<< " reaction...\n";
    Info<< " --> The forward reaction order is "
        << chem.forwardReactionOrder(r) << "\n";
    Info<< " --> The backward reaction order is "
        << chem.backwardReactionOrder(r) << "\n" << endl;

    //-------------------------------------------------------------------------
    Info<< " c-o This reaction contains a third body partner? ";
    if (chem.TBR(r)) { Info<< " YES...\n"; }
    if (!chem.TBR(r)) { Info<< " NO...\n"; }
    Info<< endl;

    //- Reaction number 
    r = 16;

    //-------------------------------------------------------------------------
    Info<< " c-o Checking reaction " << chem.elementarReaction(r) << "\n";
    Info<< " --> The forward reaction order is "
        << chem.forwardReactionOrder(r) << "\n";
    Info<< " --> The backward reaction order is "
        << chem.backwardReactionOrder(r) << "\n" << endl;

    r = 5;

    //-------------------------------------------------------------------------
    Info<< " c-o Checking reaction " << chem.elementarReaction(r) << "\n";
    Info<< " --> This reaction contains a third body partner? ";
    if (chem.TBR(r)) { Info<< " YES...\n"; }
    if (!chem.TBR(r)) { Info<< " NO...\n"; }
    Info<< " --> The collision partner is: " << chem.collisionPartner(r) <<"\n";
    Info<< " --> The reaction has a Lindemann fall off (LOW)? ";
    if (chem.LOW(r)) { Info<< " YES...\n"; }
    if (!chem.LOW(r)) { Info<< " NO...\n"; }
    Info<< " --> The fall off coefficients for the arrhenius equations are:\n --> ";
    forAll(chem.LOWCoeffs(r), coeffs)
    {
        Info<< coeffs << "   ";
    }
    Info<< endl;
    Info<< " --> The reaction does contain enhancements for the coll. partner? ";
    if (chem.ENHANCED(r)) { Info<< " YES...\n"; }
    if (!chem.ENHANCED(r)) { Info<< " NO...\n"; }

    if (chem.ENHANCED(r))
    {
        Info<< " --> For the collision partner one takes into account the "
            << "following species modifications:\n";

        const map<word, scalar>& collPartner = chem.ENHANCEDCoeffs(r);

        loopMap(species, coeff, collPartner)
        {
            Info<< "    " << species << " --> " << coeff << "\n";
        }
    }
    Info<< endl;
    Info<< " --> The reaction uses a TROE extension? ";
    if (chem.TROE(r)) { Info<< " YES...\n"; }
    if (!chem.TROE(r)) { Info<< " NO...\n"; }
    Info<< "\n" << endl;

    r = 4;

    //-------------------------------------------------------------------------
    Info<< " c-o Checking reaction " << chem.elementarReaction(r) << "\n";
    Info<< " --> This reaction contains a third body partner? ";
    if (chem.TBR(r)) { Info<< " YES...\n"; }
    if (!chem.TBR(r)) { Info<< " NO...\n"; }
    Info<< " --> The collision partner is: " << chem.collisionPartner(r) <<"\n";
    Info<< " --> The reaction has a Lindemann fall off (LOW)? ";
    if (chem.LOW(r)) { Info<< " YES...\n"; }
    if (!chem.LOW(r)) { Info<< " NO...\n"; }
    Info<< " --> The fall off coefficients for the arrhenius equations are:\n --> ";
    forAll(chem.LOWCoeffs(r), coeffs)
    {
        Info<< coeffs << "   ";
    }
    Info<< endl;
    Info<< " --> The reaction does contain enhancements for the coll. partner? ";
    if (chem.ENHANCED(r)) { Info<< " YES...\n"; }
    if (!chem.ENHANCED(r)) { Info<< " NO...\n"; }

    if (chem.ENHANCED(r))
    {
        Info<< " --> For the collision partner one takes into account the "
            << "following species modifications:\n";

        const map<word, scalar>& collPartner = chem.ENHANCEDCoeffs(r);

        loopMap(species, coeff, collPartner)
        {
            Info<< "    " << species << " --> " << coeff << "\n";
        }
    }
    Info<< endl;
    Info<< " --> The reaction uses a TROE extension? ";
    if (chem.TROE(r)) { Info<< " YES...\n"; }
    if (!chem.TROE(r)) { Info<< " NO...\n"; }
    Info<< "\n" << endl;

    r = 22;

    //-------------------------------------------------------------------------
    Info<< " c-o Checking reaction " << chem.elementarReaction(r) << "\n";
    Info<< " --> This reaction contains a third body partner? ";
    if (chem.TBR(r)) { Info<< " YES...\n"; }
    if (!chem.TBR(r)) { Info<< " NO...\n"; }
    Info<< " --> The collision partner is: " << chem.collisionPartner(r) <<"\n";
    Info<< " --> The reaction has a Lindemann fall off (LOW)? ";
    if (chem.LOW(r)) { Info<< " YES...\n"; }
    if (!chem.LOW(r)) { Info<< " NO...\n"; }
    Info<< " --> The fall off coefficients for the arrhenius equations are:\n --> ";
    forAll(chem.LOWCoeffs(r), coeffs)
    {
        Info<< coeffs << "   ";
    }
    Info<< endl;
    Info<< " --> The reaction does contain enhancements for the coll. partner? ";
    if (chem.ENHANCED(r)) { Info<< " YES...\n"; }
    if (!chem.ENHANCED(r)) { Info<< " NO...\n"; }

    if (chem.ENHANCED(r))
    {
        Info<< " --> For the collision partner one takes into account the "
            << "following species modifications:\n";

        const map<word, scalar>& collPartner = chem.ENHANCEDCoeffs(r);

        loopMap(species, coeff, collPartner)
        {
            Info<< "    " << species << " --> " << coeff << "\n";
        }
    }
    Info<< endl;
    Info<< " --> The reaction uses a TROE extension? ";
    if (chem.TROE(r)) { Info<< " YES...\n"; }
    if (!chem.TROE(r)) { Info<< " NO...\n"; }
    Info<< " --> The TROE coefficients T**, T*, T*** are:\n --> ";
    forAll(chem.TROECoeffs(r), coeffs)
    {
        Info<< coeffs << "   ";
    }
    Info<< "\n" << endl;


    //-------------------------------------------------------------------------
    for (r = 0; r < 5; ++r)
    {
        Info<< " c-o For " << chem.elementarReaction(r) << " and 1400 K, the\n"
            << " --> forward reaction rate is " << chem.kf(r, 1400) << "\n"
            << " --> backward reaction rate is " << chem.kb(r, 1400) << "\n"
            << endl;
    }


    Info<< endl;

    Footer(startTime);

    return 0;
}


// ************************************************************************* //
