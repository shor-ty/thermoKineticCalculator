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

#include "chemistryData.hpp"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

AFC::ChemistryData::ChemistryData()
:
    thermo_{false}
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

AFC::ChemistryData::~ChemistryData()
{}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //

void AFC::ChemistryData::setThermo()
{
    thermo_ = true;
}


bool AFC::ChemistryData::thermo()
{
    return thermo_;
}

// * * * * * * * * * Insert functions from ChemistryReader:: * * * * * * * * //

void AFC::ChemistryData::elements(const word element)
{
    elements_.push_back(element);
}


void AFC::ChemistryData::species(const word species)
{
    species_.push_back(species);

    //- Increment matrix
    reactionI_[species] = List<int>(0);

    //- Increment omega field
    omega_.push_back(scalar(0));
}


void AFC::ChemistryData::educt(const word species)
{
    educts_[educts_.size()-1].push_back(species);
}


void AFC::ChemistryData::product(const word species)
{
    products_[products_.size()-1].push_back(species);
}


void AFC::ChemistryData::nuEducts(const word species, const int nu)
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


void AFC::ChemistryData::nuProducts(const word species, const int nu)
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


void AFC::ChemistryData::elementarReaction(const string reaction)
{
    elementarReaction_[nReac_] = reaction;
}


void AFC::ChemistryData::ignoredElementarReaction(const string reaction)
{
    ignoredElementarReaction_[nIgnored_-1] = reaction;
}


void AFC::ChemistryData::arrheniusCoeffs
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


void AFC::ChemistryData::LOWCoeffs
(
    const scalar coeff,
    const unsigned int coeffNo
)
{

    LOWCoeffs_[nReac_][coeffNo] = coeff;
}


void AFC::ChemistryData::TROECoeffs
(
    const scalar coeff,
    const unsigned int coeffNo
)
{
    TROECoeffs_[nReac_][coeffNo] = coeff;
}


void AFC::ChemistryData::SRICoeffs
(
    const scalar coeff,
    const unsigned int coeffNo
)
{
    SRICoeffs_[nReac_][coeffNo] = coeff;
}


void AFC::ChemistryData::ENHANCEDCoeffs
(
    const word species,
    const scalar factor
)
{
    ENHANCEDCoeffs_[nReac_][species] = factor;
}


void AFC::ChemistryData::incrementDuplicated()
{
    nDuplicated_++;
}


void AFC::ChemistryData::incrementIgnored()
{
    nIgnored_++;

    //- In this regard, we also increment the container
    ignoredElementarReaction_.push_back("");
}


void AFC::ChemistryData::incrementReac()
{
    nReac_++;
}


void AFC::ChemistryData::speciesInReaction(const word species)
{
    speciesInReaction_[nReac_].push_back(species);
}


void AFC::ChemistryData::incrementMatrixesVectors()
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

    //- boolList for backward reaction
    backwardReaction_.push_back(false);

    //- List<List<word> >
    speciesInReaction_.push_back(List<word>(0));

    //- Matrix of Arrhenius coeffs
    arrheniusCoeffs_.push_back(scalarField(3, scalar(0)));

    //- Matrix of ARRHENIUS coeffs for LOW pressure
    LOWCoeffs_.push_back(scalarField(3, scalar(0)));

    //- Matrix of TROE coeffs
    TROECoeffs_.push_back(scalarField(4, scalar(0)));

    //- Matrix of SRI coeffs
    SRICoeffs_.push_back(scalarField(5));

    //- MapList of enhanced factors for adjustment (species + value)
    ENHANCEDCoeffs_.push_back(map<word,scalar>());

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

void AFC::ChemistryData::BR(const bool set)
{
    backwardReaction_[nReac_] = set;
}


void AFC::ChemistryData::TBR(const bool set)
{
    TBR_[nReac_] = set;
}


void AFC::ChemistryData::LOW(const bool set)
{
    LOW_[nReac_] = set;
}


void AFC::ChemistryData::TROE(const bool set)
{
    TROE_[nReac_] = set;
}


void AFC::ChemistryData::SRI(const bool set)
{
    SRI_[nReac_] = set;
}


void AFC::ChemistryData::ENHANCE(const bool set)
{
    ENHANCE_[nReac_] = set;
}


void AFC::ChemistryData::setReacNumbers(const word species, const int r)
{
    reactionI_.at(species).push_back(r);
}


void AFC::ChemistryData::forwardReactionOrder()
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


void AFC::ChemistryData::backwardReactionOrder()
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

void AFC::ChemistryData::updateKf(const int r, const scalar kf)
{
    kf_[r] = kf;
}


void AFC::ChemistryData::updateKb(const int r, const scalar kb)
{
    kb_[r] = kb;
}


void AFC::ChemistryData::updateKc(const int r, const scalar Kc)
{
    Kc_[r] = Kc;
}


void AFC::ChemistryData::calculateOmega(const int s, const scalar omega)
{
    omega_[s] = omega;
}


void AFC::ChemistryData::calculateOmega(const scalarField& omega)
{
    omega_ = omega;
}


void AFC::ChemistryData::updateGlobalReactionOrder()
{
    //- Forward reaction order
    const scalar& fRO = forwardReactionOrder(nReac_);

    //- Backward reaction order
    const scalar& bRO = backwardReactionOrder(nReac_);

    //- Global reaction order
    globalReactionOrder_[nReac_] = fRO + bRO;
}

// * * * * * * * * * * * * * * * Return functions  * * * * * * * * * * * * * //

bool AFC::ChemistryData::BR(const int reacNo) const
{
    return backwardReaction_[reacNo];
}


bool AFC::ChemistryData::BR(const int reacNo)
{
    return backwardReaction_[reacNo];
}


bool AFC::ChemistryData::TBR(const int reacNo) const
{
    return TBR_[reacNo];
}


bool AFC::ChemistryData::LOW(const int reacNo) const
{
    return LOW_[reacNo];
}


bool AFC::ChemistryData::TROE(const int reacNo) const
{
    return TROE_[reacNo];
}


bool AFC::ChemistryData::SRI(const int reacNo) const
{
    return SRI_[reacNo];
}


bool AFC::ChemistryData::ENHANCED(const int reacNo) const
{
    return ENHANCE_[reacNo];
}


AFC::scalar AFC::ChemistryData::dS() const
{
    return dS_;
}


AFC::scalar AFC::ChemistryData::dG() const
{
    return dG_;
}


AFC::wordList AFC::ChemistryData::elements() const
{
    return elements_;
}


AFC::wordList AFC::ChemistryData::species() const
{
    return species_;
}


AFC::wordList AFC::ChemistryData::speciesEducts(const int r) const
{
    return educts_[r];
}


AFC::wordList AFC::ChemistryData::speciesProducts(const int r) const
{
    return products_[r];
}


AFC::map<AFC::word, int> AFC::ChemistryData::nuEducts(const int r ) const
{
    return nuEducts_[r];
}


AFC::map<AFC::word, int> AFC::ChemistryData::nuProducts(const int r) const
{
    return nuProducts_[r];
}


unsigned int AFC::ChemistryData::nDublicated() const
{
    return nDuplicated_;
}


unsigned int AFC::ChemistryData::nIgnored() const
{
    return nIgnored_;
}


int AFC::ChemistryData::nReac() const
{
    return nReac_;
}


AFC::List<AFC::string> AFC::ChemistryData::elementarReaction() const
{
    return elementarReaction_;
}


AFC::string AFC::ChemistryData::elementarReaction(const int r) const
{
    return elementarReaction_[r];
}


AFC::List<AFC::string> AFC::ChemistryData::ignoredElementarReaction() const
{
    return ignoredElementarReaction_;
}


AFC::string AFC::ChemistryData::ignoredElementarReaction(const int r) const
{
    return ignoredElementarReaction_[r];
}


AFC::List<int> AFC::ChemistryData::reacNumbers(const word species) const
{
    return reactionI_.at(species);
}


AFC::wordMatrix AFC::ChemistryData::speciesInReaction() const
{
    return speciesInReaction_;
}


AFC::wordList AFC::ChemistryData::speciesInReaction(const int r) const
{
    return speciesInReaction_[r];
}


AFC::scalarList AFC::ChemistryData::arrheniusCoeffs(const int reacNo) const
{
    return arrheniusCoeffs_[reacNo];
}


AFC::scalarList AFC::ChemistryData::LOWCoeffs(const int reacNo) const
{
    return LOWCoeffs_[reacNo];
}


AFC::scalarList AFC::ChemistryData::TROECoeffs(const int reacNo) const
{
    return TROECoeffs_[reacNo];
}


AFC::scalarList AFC::ChemistryData::SRICoeffs(const int reacNo) const
{
    return SRICoeffs_[reacNo];
}


AFC::map<AFC::word, AFC::scalar>
AFC::ChemistryData::ENHANCEDCoeffs(const int reacNo) const
{
    return ENHANCEDCoeffs_[reacNo];
}


AFC::scalar AFC::ChemistryData::kf(const int reacNo) const
{
    return kf_[reacNo];
}


AFC::scalarList AFC::ChemistryData::kf() const
{
    return kf_;
}


AFC::scalar AFC::ChemistryData::kb(const int reacNo) const
{
    return kb_[reacNo];
}


AFC::scalarList AFC::ChemistryData::kb() const
{
    return kb_;
}


AFC::scalar AFC::ChemistryData::Kc(const int reacNo) const
{
    return Kc_[reacNo];
}


AFC::scalarList AFC::ChemistryData::Kc() const
{
    return Kc_;
}


AFC::scalar AFC::ChemistryData::omega(const int s) const
{
    return omega_[s];
}


AFC::scalarField AFC::ChemistryData::omega() const
{
    return omega_;
}


AFC::scalar AFC::ChemistryData::forwardReactionOrder(const int r) const
{
    return forwardReactionOrder_[r];
}


AFC::scalar AFC::ChemistryData::backwardReactionOrder(const int r) const
{
    return backwardReactionOrder_[r];
}


AFC::scalar AFC::ChemistryData::globalReactionOrder(const int r) const
{
    return globalReactionOrder_[r];
}


// ************************************************************************* //
