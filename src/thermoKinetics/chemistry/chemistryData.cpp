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

#include "chemistryData.hpp"
#include "chemistryReader.hpp"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

TKC::ChemistryData::ChemistryData(const string fileName)
//:
    //thermo_{false}
{
    ChemistryReader chemReader(fileName);

    chemReader.read(*this);
}



// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

TKC::ChemistryData::~ChemistryData()
{}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //

/*void TKC::ChemistryData::setThermo()
{
    thermo_ = true;
}


bool TKC::ChemistryData::thermo() const
{
    return thermo_;
}
*/


// * * * * * * * * * Insert functions from ChemistryReader:: * * * * * * * * //

void TKC::ChemistryData::elements(const word element)
{
    elements_.push_back(element);
}


void TKC::ChemistryData::species(const word species)
{
    species_.push_back(species);

    //- Increment matrix
    reactionI_[species] = List<int>(0);

    //- Increment omega field
    omega_.push_back(scalar(0));
}


void TKC::ChemistryData::educt(const word species)
{
    educts_[educts_.size()-1].push_back(species);
}


void TKC::ChemistryData::product(const word species)
{
    products_[products_.size()-1].push_back(species);
}


void TKC::ChemistryData::nuEducts(const word species, const int nu)
{
    //- Check if species already there, if yes, increment nu
    if (nuEducts_[nReac_].find(species) != nuEducts_[nReac_].end())
    {
        nuEducts_[nReac_][species] += nu;
    }
    else
    {
        nuEducts_[nReac_][species] = nu;
    }
}


void TKC::ChemistryData::nuProducts(const word species, const int nu)
{
    //- Check if species already there, if yes, increment nu
    if (nuProducts_[nReac_].find(species) != nuProducts_[nReac_].end())
    {
        nuProducts_[nReac_][species] += nu;
    }
    else
    {
        nuProducts_[nReac_][species] = nu;
    }
}


void TKC::ChemistryData::elementarReaction(const string reaction)
{
    elementarReaction_[nReac_] = reaction;
}


void TKC::ChemistryData::duplicatedElementarReaction(const string reaction)
{
    duplicatedElementarReaction_[nReac_] = reaction;
}


void TKC::ChemistryData::ignoredElementarReaction(const string reaction)
{
    ignoredElementarReaction_[nIgnored_-1] = reaction;
}


void TKC::ChemistryData::arrheniusCoeffs
(
    const scalar coeff_1,
    const scalar coeff_2,
    const scalar coeff_3
)
{
    arrheniusCoeffs_[nReac_][0] = coeff_1;
    arrheniusCoeffs_[nReac_][1] = coeff_2;
    arrheniusCoeffs_[nReac_][2] = coeff_3;
}


void TKC::ChemistryData::collisionPartner(const word collisionPartner)
{
    collisionPartner_[nReac_] = collisionPartner;
}


void TKC::ChemistryData::LOWCoeffs
(
    const scalar coeff,
    const unsigned int coeffNo
)
{

    LOWCoeffs_[nReac_][coeffNo] = coeff;
}


void TKC::ChemistryData::TROECoeffs
(
    const scalar coeff,
    const unsigned int coeffNo
)
{
    TROECoeffs_[nReac_][coeffNo] = coeff;
}


void TKC::ChemistryData::SRICoeffs
(
    const scalar coeff,
    const unsigned int coeffNo
)
{
    SRICoeffs_[nReac_][coeffNo] = coeff;
}


void TKC::ChemistryData::ENHANCEDCoeffs
(
    const word species,
    const scalar factor
)
{
    ENHANCEDCoeffs_[nReac_].at(species) = factor;
}


void TKC::ChemistryData::incrementDuplicated()
{
    nDuplicated_++;
    duplicatedElementarReaction_.push_back("");
}


void TKC::ChemistryData::incrementIgnored()
{
    nIgnored_++;
    ignoredElementarReaction_.push_back("");
}


void TKC::ChemistryData::incrementReac()
{
    nReac_++;
}


void TKC::ChemistryData::speciesInReaction(const word species)
{
    speciesInReaction_[nReac_].push_back(species);
}


void TKC::ChemistryData::incrementMatrixesVectors()
{
    //- stringList for saving reactions
    elementarReaction_.push_back("");

    //- List<List<word> > for educts
    educts_.push_back(List<word>(0));

    //- List<List<word> > for products
    products_.push_back(List<word>(0));

    //- MapList for for stochiometric factors for educts (species + value)
    nuEducts_.push_back(map<word, int>());

    //- MapList for for stochiometric factors for products (species + value)
    nuProducts_.push_back(map<word, int>());

    //- scalarList for forward reaction order
    forwardReactionOrder_.push_back(0);

    //- scalarList for backward reaction order
    backwardReactionOrder_.push_back(0);

    //- scalarList for global reaction order
    globalReactionOrder_.push_back(0);

    //- Collision partner
    collisionPartner_.push_back("");

    //- boolList for THIRD BODY REACTION
    TBR_.push_back(false);

    //- boolList for THIRD BODY REACTION LOW
    LOW_.push_back(false);

    //- boolList for THIRD BODY REACTION TROE
    TROE_.push_back(false);

    //- boolList for THIRD BODY REACTION SRI
    SRI_.push_back(false);

    //- boolList for THIRD BODY REACTION of ENHANCEMENT FACTORS
    ENHANCE_.push_back(false);

    //- boolList for forward reaction
    forwardReaction_.push_back(false);

    //- boolList for backward reaction
    backwardReaction_.push_back(false);

    //- List<List<word> >
    speciesInReaction_.push_back(List<word>(0));

    //- List of Arrhenius coeffs
    arrheniusCoeffs_.push_back(scalarField(3, scalar(0)));

    //- List of Arrhenius coeffs for Low pressure (fall off)
    LOWCoeffs_.push_back(scalarField(3, scalar(0)));

    //- List of TROE coeffs
    TROECoeffs_.push_back(scalarField(4, scalar(0)));

    //- List of SRI coeffs
    SRICoeffs_.push_back(scalarField(5));

    //- MapList of enhanced factors for adjustment (species + value)
    //  For each reaction we save it
    {
        map<word, scalar> tmp;

        forAll(species(), s)
        {
            tmp[s] = 1;
        }

        //- Add all species with factor = 1
        ENHANCEDCoeffs_.push_back(tmp);
    }

    //- Update entrys in SRI, d and e if not used
    {
        SRICoeffs_[nReac_][3] = 1;
        SRICoeffs_[nReac_][4] = 0;
    }

    //- Reaction rates kf
    kf_.push_back(scalar(0));

    //- Reaction rates kb
    kb_.push_back(scalar(0));

    //- Equilibrium constant Kc
    Kc_.push_back(scalar(0));
}


// * * * * * * * * * * * * * Setter bool functions * * * * * * * * * * * * * //

void TKC::ChemistryData::FR(const bool set)
{
    forwardReaction_[nReac_] = set;
}


void TKC::ChemistryData::BR(const bool set)
{
    backwardReaction_[nReac_] = set;
}


void TKC::ChemistryData::TBR(const bool set)
{
    TBR_[nReac_] = set;
}


void TKC::ChemistryData::LOW(const bool set)
{
    LOW_[nReac_] = set;
}


void TKC::ChemistryData::TROE(const bool set)
{
    TROE_[nReac_] = set;
}


void TKC::ChemistryData::SRI(const bool set)
{
    SRI_[nReac_] = set;
}


void TKC::ChemistryData::ENHANCE(const bool set)
{
    ENHANCE_[nReac_] = set;
}


void TKC::ChemistryData::setReacNumbers(const word species, const int r)
{
    reactionI_.at(species).push_back(r);
}


void TKC::ChemistryData::forwardReactionOrder()
{
    //- Stochiometric coefficients of educt species
    const map<word, int>& educt = nuEducts(nReac_);

    scalar order{0};

    forAll(educt, e)
    {
        order += e.second;
    }

    //- Set the reaction order
    forwardReactionOrder_[nReac_] = order;
}


void TKC::ChemistryData::backwardReactionOrder()
{
    //- Stochiometric coefficients of product species
    const map<word, int>& product = nuProducts(nReac_);

    scalar order{0};

    forAll(product, p)
    {
        order += p.second;
    }

    //- Set the reaction order
    backwardReactionOrder_[nReac_] = order;
}

// * * * * * * * * * * * * * * * Update functions  * * * * * * * * * * * * * //

void TKC::ChemistryData::updateKf(const int r, const scalar kf)
{
    kf_[r] = kf;
}


void TKC::ChemistryData::updateKb(const int r, const scalar kb)
{
    kb_[r] = kb;
}


void TKC::ChemistryData::updateKc(const int r, const scalar Kc)
{
    Kc_[r] = Kc;
}


void TKC::ChemistryData::calculateOmega(const int s, const scalar omega)
{
    omega_[s] = omega;
}


void TKC::ChemistryData::calculateOmega(const scalarField& omega)
{
    omega_ = omega;
}


void TKC::ChemistryData::updateGlobalReactionOrder()
{
    //- Forward reaction order
    const scalar& fRO = forwardReactionOrder(nReac_);

    //- Backward reaction order
    const scalar& bRO = backwardReactionOrder(nReac_);

    //- Global reaction order
    globalReactionOrder_[nReac_] = fRO + bRO;
}


// * * * * * * * * * * * * * * * Return functions  * * * * * * * * * * * * * //

bool TKC::ChemistryData::FR(const int reacNo) const
{
    return forwardReaction_[reacNo];
}


bool TKC::ChemistryData::BR(const int reacNo) const
{
    return backwardReaction_[reacNo];
}


bool TKC::ChemistryData::TBR(const int reacNo)
{
    return TBR_[reacNo];
}


bool TKC::ChemistryData::TBR(const int reacNo) const
{
    return TBR_[reacNo];
}


bool TKC::ChemistryData::LOW(const int reacNo) const
{
    return LOW_[reacNo];
}


bool TKC::ChemistryData::TROE(const int reacNo) const
{
    return TROE_[reacNo];
}


bool TKC::ChemistryData::SRI(const int reacNo) const
{
    return SRI_[reacNo];
}


bool TKC::ChemistryData::ENHANCED(const int reacNo) const
{
    return ENHANCE_[reacNo];
}


TKC::scalar TKC::ChemistryData::dS() const
{
    return dS_;
}


TKC::scalar TKC::ChemistryData::dG() const
{
    return dG_;
}


TKC::wordList TKC::ChemistryData::elements() const
{
    return elements_;
}


TKC::wordList TKC::ChemistryData::species() const
{
    return species_;
}


TKC::wordList TKC::ChemistryData::educts(const int r) const
{
    return educts_[r];
}



TKC::wordList TKC::ChemistryData::products(const int r) const
{
    return products_[r];
}


TKC::map<TKC::word, int> TKC::ChemistryData::nuEducts(const int r) const
{
    return nuEducts_[r];
}


TKC::map<TKC::word, int> TKC::ChemistryData::nuProducts(const int r) const
{
    return nuProducts_[r];
}


unsigned int TKC::ChemistryData::nDuplicated() const
{
    return nDuplicated_;
}


unsigned int TKC::ChemistryData::nIgnored() const
{
    return nIgnored_;
}


int TKC::ChemistryData::nReac() const
{
    return nReac_;
}


TKC::stringList TKC::ChemistryData::elementarReaction() const
{
    return elementarReaction_;
}


TKC::string TKC::ChemistryData::elementarReaction(const int r) const
{
    return elementarReaction_[r];
}


TKC::List<TKC::string> TKC::ChemistryData::ignoredElementarReaction() const
{
    return ignoredElementarReaction_;
}


TKC::string TKC::ChemistryData::ignoredElementarReaction(const int r) const
{
    return ignoredElementarReaction_[r];
}


TKC::List<int> TKC::ChemistryData::reacNumbers(const word species) const
{
    return reactionI_.at(species);
}


TKC::List<TKC::wordList> TKC::ChemistryData::speciesInReaction() const
{
    return speciesInReaction_;
}


TKC::wordList TKC::ChemistryData::speciesInReaction(const int r) const
{
    return speciesInReaction_[r];
}


bool TKC::ChemistryData::forwardReaction(const int r) const
{
    return forwardReaction_[r];
}


bool TKC::ChemistryData::backwardReaction(const int r) const
{
    return backwardReaction_[r];
}


TKC::scalarList TKC::ChemistryData::arrheniusCoeffs(const int reacNo) const
{
    return arrheniusCoeffs_[reacNo];
}


TKC::word TKC::ChemistryData::collisionPartner(const int reacNo) const
{
    return collisionPartner_[reacNo];
}


TKC::scalarList TKC::ChemistryData::LOWCoeffs(const int reacNo) const
{
    return LOWCoeffs_[reacNo];
}


TKC::scalarList TKC::ChemistryData::TROECoeffs(const int reacNo) const
{
    return TROECoeffs_[reacNo];
}


TKC::scalarList TKC::ChemistryData::SRICoeffs(const int reacNo) const
{
    return SRICoeffs_[reacNo];
}


TKC::map<TKC::word, TKC::scalar>
TKC::ChemistryData::ENHANCEDCoeffs(const int reacNo) const
{
    return ENHANCEDCoeffs_[reacNo];
}


TKC::scalar TKC::ChemistryData::kf(const int reacNo) const
{
    return kf_[reacNo];
}


TKC::scalarList TKC::ChemistryData::kf() const
{
    return kf_;
}


TKC::scalar TKC::ChemistryData::kb(const int reacNo) const
{
    return kb_[reacNo];
}


TKC::scalarList TKC::ChemistryData::kb() const
{
    return kb_;
}


TKC::scalar TKC::ChemistryData::Kc(const int reacNo) const
{
    return Kc_[reacNo];
}


TKC::scalarList TKC::ChemistryData::Kc() const
{
    return Kc_;
}


TKC::scalar TKC::ChemistryData::omega(const int s) const
{
    return omega_[s];
}


TKC::scalarField TKC::ChemistryData::omega() const
{
    return omega_;
}


TKC::scalar TKC::ChemistryData::forwardReactionOrder(const int r) const
{
    return forwardReactionOrder_[r];
}


TKC::scalar TKC::ChemistryData::backwardReactionOrder(const int r) const
{
    return backwardReactionOrder_[r];
}


TKC::scalar TKC::ChemistryData::globalReactionOrder(const int r) const
{
    return globalReactionOrder_[r];
}


// ************************************************************************* //
