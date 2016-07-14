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
    nReac_{-1},

    thermo_{false}
{
    if (debug)
    {
        Info<< "Constructor ChemistryData\n";
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

AFC::ChemistryData::~ChemistryData()
{
    if (debug)
    {
        Info<< "Destructor ChemistryData\n";
    }
}


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

void AFC::ChemistryData::elements
(
    const word& element
)
{
    elements_.push_back(element);
}


void AFC::ChemistryData::species
(
    const word& species
)
{
    species_.push_back(species);
    
    //- Increment matrix
    reactionI_[species] = List<int>(0);

    //- Increment omega field
    omega_.push_back(scalar(0));
}


void AFC::ChemistryData::elementarReaction
(
    const string& reaction
)
{
    if (debug)
    {
        Info<< " --> AFC::ChemistryData::elementarReaction" << endl;
    }

    elementarReaction_[nReac_] = reaction;
}


void AFC::ChemistryData::arrheniusCoeffs
(
    const scalar& coeff_1,
    const scalar& coeff_2,
    const scalar& coeff_3
)
{
    if (debug)
    {
        Info<< " --> AFC::ChemistryData::arrheniusCoeffs" << endl;
    }

    arrheniusCoeffs_[nReac_][0] = coeff_1;
    arrheniusCoeffs_[nReac_][1] = coeff_2;
    arrheniusCoeffs_[nReac_][2] = coeff_3;
}


void AFC::ChemistryData::LOWCoeffs
(
    const scalar& coeff,
    const unsigned int& coeffNo
)
{
    if (debug)
    {
        Info<< " --> AFC::ChemistryData::LOWCoeffs" << endl;
    }

    LOWCoeffs_[nReac_][coeffNo] = coeff;
}


void AFC::ChemistryData::TROECoeffs
(
    const scalar& coeff,
    const unsigned int& coeffNo
)
{
    if (debug)
    {
        Info<< " --> AFC::ChemistryData::TROECoeffs" << endl;
    }

    TROECoeffs_[nReac_][coeffNo] = coeff;
}


void AFC::ChemistryData::SRICoeffs
(
    const scalar& coeff,
    const unsigned int& coeffNo
)
{
    if (debug)
    {
        Info<< " --> AFC::ChemistryData::SRICoeffs" << endl;
    }

    SRICoeffs_[nReac_][coeffNo] = coeff;
}


void AFC::ChemistryData::enhanceFactors
(
    const word& species,
    const scalar& value
)
{
    if (debug)
    {
        Info<< " --> AFC::ChemistryData::enhanceFactors" << endl;
    }

    enhancedFactors_[nReac_][species] = value;
}


void AFC::ChemistryData::incrementDuplicated()
{
    nDuplicated_++;
}


void AFC::ChemistryData::incrementReac()
{
    nReac_++;
}


void AFC::ChemistryData::speciesNu
(
    const word& species,
    const scalar& nu
)
{
    nu_[nReac_][species] = nu;
}


void AFC::ChemistryData::speciesInReaction
(
    const word& species
)
{
    speciesInReaction_[nReac_].push_back(species);
}


void AFC::ChemistryData::incrementMatrixesVectors()
{
    if (debug)
    {
        Info<< " --> AFC::ChemistryData::incrementMatrixesVectors" << endl;
    }

    //- stringList for saving reactions
    elementarReaction_.push_back("");

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

    //- MapList of enhanced factors (species + value)
    enhancedFactors_.push_back(map<word,scalar>());

    //- MapList of stochiometric factors (species + value)
    nu_.push_back(map<word,scalar>());

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


void AFC::ChemistryData::setReacNumbers
(
    const word& species,
    const int& r
)
{
    reactionI_.at(species).push_back(r);
}


// * * * * * * * * * * * * * * * Update functions  * * * * * * * * * * * * * //

void AFC::ChemistryData::updateKf
(
    const int& r,
    const scalar& kf
)
{
    kf_[r] = kf;
}


void AFC::ChemistryData::updateKb
(
    const int& r,
    const scalar& kb
)
{
    kb_[r] = kb;
}


void AFC::ChemistryData::updateKc
(
    const int& r,
    const scalar& Kc
)
{
    Kc_[r] = Kc;
}


void AFC::ChemistryData::calculateOmega
(
    const int& s,
    const scalar& omega
)
{
    omega_[s] = omega;
}


void AFC::ChemistryData::calculateOmega
(
    const scalarField& omega
)
{
    omega_ = omega;
}


// * * * * * * * * * * * * * * * Return functions  * * * * * * * * * * * * * //

bool AFC::ChemistryData::BR
(
    const int& reacNo
) const
{
    return backwardReaction_[reacNo];
}


bool AFC::ChemistryData::TBR
(
    const int& reacNo
) const
{
    return TBR_[reacNo];
}


bool AFC::ChemistryData::LOW
(
    const int& reacNo
) const
{
    return LOW_[reacNo];
}


bool AFC::ChemistryData::TROE
(
    const int& reacNo
) const
{
    return TROE_[reacNo];
}


bool AFC::ChemistryData::SRI
(
    const int& reacNo
) const
{
    return SRI_[reacNo];
}


bool AFC::ChemistryData::ENHANCED
(
    const int& reacNo
) const
{
    return ENHANCE_[reacNo];
}


AFC::scalar AFC::ChemistryData::dH() const
{
    return dH_;
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


int AFC::ChemistryData::nReac() const
{
    return nReac_;
}


AFC::List<AFC::string> AFC::ChemistryData::elementarReaction() const
{
    return elementarReaction_;
}


AFC::string AFC::ChemistryData::elementarReaction
(
    const int& r
) const
{
    return elementarReaction_[r];
}


AFC::List<int> AFC::ChemistryData::reacNumbers
(
    const word& species
) const
{
    return reactionI_.at(species);
}


AFC::wordMatrix AFC::ChemistryData::speciesInReaction() const
{
    return speciesInReaction_;
}


AFC::wordList AFC::ChemistryData::speciesInReaction
(
    const int& r 
) const
{
    return speciesInReaction_[r];
}


AFC::scalarList AFC::ChemistryData::arrheniusCoeffs
(
    const int& reacNo
) const
{
    return arrheniusCoeffs_[reacNo];
}


AFC::scalarList AFC::ChemistryData::LOWCoeffs
(
    const int& reacNo
) const
{
    return LOWCoeffs_[reacNo];
}


AFC::scalarList AFC::ChemistryData::TROECoeffs
(
    const int& reacNo
) const
{
    return TROECoeffs_[reacNo];
}


AFC::scalarList AFC::ChemistryData::SRICoeffs
(
    const int& reacNo
) const
{
    return SRICoeffs_[reacNo];
}


AFC::wordList AFC::ChemistryData::enhancedSpecies
(
    const int& reacNo
) const
{
    wordList enhancedSpecies;

    //- Need map iterator
    map<word, scalar>::iterator iter;
    map<word, scalar> tmp = enhancedFactors_[reacNo];

    for
    (
        iter = tmp.begin();
        iter != tmp.end();
        iter++
    )
    {
        enhancedSpecies.push_back((*iter).first);
    }

    return enhancedSpecies;
}


AFC::map<AFC::word, AFC::scalar> AFC::ChemistryData::enhancedFactors
(
    const int& reacNo
) const
{
    return enhancedFactors_[reacNo];
}


AFC::scalar AFC::ChemistryData::enhancedFactors
(
    const int& reacNo,
    const word& species
) const
{
    return enhancedFactors_[reacNo].at(species);
}


AFC::scalar AFC::ChemistryData::kf
(
    const int& reacNo 
) const
{
    return kf_[reacNo];
}


AFC::scalarList AFC::ChemistryData::kf() const
{
    return kf_;
}


AFC::scalar AFC::ChemistryData::kb
(
    const int& reacNo 
) const
{
    return kb_[reacNo];
}


AFC::scalarList AFC::ChemistryData::kb() const
{
    return kb_;
}


AFC::scalar AFC::ChemistryData::Kc
(
    const int& reacNo
) const
{
    return Kc_[reacNo];
}


AFC::scalarList AFC::ChemistryData::Kc() const
{
    return Kc_;
}


AFC::scalar AFC::ChemistryData::omega
(
    const int& s
) const
{
    return omega_[s];
}


AFC::scalarField AFC::ChemistryData::omega() const
{
    return omega_;
}


AFC::map<AFC::word, AFC::scalar> AFC::ChemistryData::nu
(
    const int& r 
) const
{
    return nu_[r];
}


AFC::mapList<AFC::word, AFC::scalar> AFC::ChemistryData::nu() const
{
    return nu_;
}


// ************************************************************************* //
