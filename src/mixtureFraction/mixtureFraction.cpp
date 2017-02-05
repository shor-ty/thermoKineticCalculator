/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2015 Holzmann-cfd
  c     c-o-o-o             |
-------------------------------------------------------------------------------
License This file is part of Automatic Flamelet Constructor.

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

#include "mixtureFraction.hpp"
#include "constants.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::MixtureFraction::MixtureFraction
(
    Chemistry& chem,
    const Thermo& therm,
    const Transport& trans,
    const Properties& prop,
    const scalar& sDR,
    const scalar& defect
)
: 
    defect_(defect),

    sDR_(sDR),

    thermo_(therm),

    transport_(trans),

    chemistry_(chem),

    properties_(prop)
{
    if (debug_)
    {
        Info<< "Constructor MixtureFraction\n" << endl;
    }

    //- Number of points for the discretisation
    //  The number also contains the boundary conditions
    const int nZ = prop.nZPoints();

    //- Calculate deltaZ
    const scalar delta = 1. / (nZ-1);

    //- Temperature of fuel and oxidizer
    const scalar& TF = prop.fuelTemperature();
    const scalar& TO = prop.oxidizerTemperature();

    //- Take the lower temperature for initialisation
    const scalar Tmin = min(TF, TO);

    //- Resize fields and init with zero
    {
        Z_.resize(nZ, 0);
        T_.resize(nZ, Tmin);
        rho_.resize(nZ, 0);
        cp_.resize(nZ, 0);
        MMW_.resize(nZ, 0);
        mu_.resize(nZ, 0);
        lambda_.resize(nZ, 0);
        Cmix_.resize(nZ, 0);
    }

    //- Resize map fields
    {
        //- Species that are used in the combustion
        const wordList& species_ = species();

        //- Generate tmp map with all species and zero value
        map<word, scalar> tmp;

        forAll(species_, s)
        {
            tmp[s] = 0; 
        }

        //- Resize and init
        {
            speciesMol_.resize(nZ, tmp);
            speciesMass_.resize(nZ, tmp);
            speciesCon_.resize(nZ, tmp);
            omega_.resize(nZ, tmp);
        }
    }


    //- Calculate values for Z and save
    for (int i=0; i<nZ; ++i)
    {
        Z_[i] = i * delta;
    }

    //- Set temperature boundary conditions
    T_[0] = prop.oxidizerTemperature();
    T_[nZ-1] = prop.fuelTemperature();


    // Initialize with mole fraction
    if (prop.input() == "mol")
    {
        //- Calculate mole fraction (initial linear distribution)
        /*{
            const wordList& species = chem.species();

            //- Set the oxidizer mol fraction
            map<word, scalar> oxidizerMolFraction = prop.oxidizerCompMol();
            //
            //- Set the fuel mol fraction
            map<word, scalar> fuelMolFraction = prop.fuelCompMol();

            //- Set all species to zero
            forAll(species, s)
            {
                if
                (
                    oxidizerMolFraction[s] > 1e-10
                 || fuelMolFraction[s] > 1e-10
                )
                {
                    speciesMol_[s] =
                        oxidizerMolFraction[s] * (1 - Zvalue)
                      + fuelMolFraction[s] * Zvalue;
                }
                else
                {
                    speciesMol_[s] = 0;
                }
            }

            if (debug)
            {
                Info<< "Discrete point Z (mol): " << Z_ << endl;

                forAll(species, s)
                {
                    if (speciesMol_.at(s) > 0)
                    {
                        Info<< s << ": " << speciesMol_.at(s) << endl;
                    }
                }
            }

            //- Calculate mass fraction Y
            XtoY();
        }*/
    }

    // Initialize with mass fraction
    else if (prop.input() == "mass")
    {
        //- Actual field is mass fraction
        actualField_ = "mass";

        //- Set boundary conditions for oxidizer stream
        {
            const wordList& speciesO = prop.speciesOxidizer();
            const map<word, scalar>& Y = prop.oxidizerY();

            forAll(speciesO, s)
            {
                speciesMass_[0][s] = Y.at(s); 
            }
        }
        //- Set boundary conditions for fuel stream 
        {
            const wordList& speciesF = prop.speciesFuel();
            const map<word, scalar>& Y = prop.fuelY();

            forAll(speciesF, s)
            {
                speciesMass_[nZ-1][s] = Y.at(s); 
            }
        }

        //- Update fields
        updateFields();

        //- Update denisty
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::MixtureFraction::~MixtureFraction()
{
    if (debug_)
    {
//        Info<< "Destruct MixtureFraction\n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void AFC::MixtureFraction::updateMeanMW()
{
    if (updatedMeanMW())
    {
        return;
    }

    const wordList& species_ = species();        

    //- Mixture fraction space
    for (int Z = 0; Z < nZPoints(); ++Z)
    {
        //- Get access
        scalar& value = MMW_[Z];

        //- Reset
        value = 0;

        //- Calculate mean moleculare weight with mol fraction
        if (actualField_ == "mol")
        {
            forAll(species_, s)
            {
                value += speciesMol_[Z].at(s) * thermo_.MW(s);
            }
        } 
        
        //- Calculate mean moleculare weight with mass fraction
        else if (actualField_ == "mass")
        {
            scalar tmp{0};

            forAll(species_, s)
            {
                tmp += speciesMass_[Z].at(s) / thermo_.MW(s);
            }

            //- If tmp == 0, no species available MW = 0
            if (tmp < 1e-6)
            {
                value = 0;
            }
            else
            {
                value = 1/tmp;
            }
        }
        //- Calculate mean moleculare weight with concentration 
        else if (actualField_ == "con")
        {
            scalar numerator{0};
            scalar denominator{0};

            forAll(species_, s)
            {
                numerator += speciesCon_[Z].at(s)*thermo_.MW(s);
                denominator += speciesCon_[Z].at(s);
            }

            value = numerator / denominator;
        }
        else
        {
            FatalError
            (
                "    Calculation of mean moleculare weight not implemented.\n"
                "    You can either use 'mass', 'mol' or 'con' as argument.",
                __FILE__,
                __LINE__
            );
        }
    }

    //- Set update switch to true
    updatedMeanMW(true); 
}


/*void AFC::MixtureFraction::calculateMeanCp
(
    const scalar& T 
)
{
    const wordList& species = chemistry_.species();        

    //- Reset
    cp_ = 0;

    forAll(species, s)
    {
        cp_ += speciesMass_.at(species[s]) * thermo_.cp(species[s], T);
    }
}


void AFC::MixtureFraction::calculateMeanH
(
    const scalar& T 
)
{
    const wordList& species = chemistry_.species();        

    //- Reset
    H_ = 0;

    forAll(species, s)
    {
        H_ += speciesMass_.at(species[s]) * thermo_.H(species[s], T);
    }
}


void AFC::MixtureFraction::calculateMeanS
(
    const scalar& T 
)
{
    const wordList& species = chemistry_.species();        

    //- Reset
    S_ = 0;

    forAll(species, s)
    {
        S_ += speciesMass_.at(species[s]) * thermo_.S(species[s], T);
    }
}


void AFC::MixtureFraction::calculateMeanG
(
    const scalar& T 
)
{
    const wordList& species = chemistry_.species();        

    //- Reset
    G_ = 0;

    forAll(species, s)
    {
        G_ += speciesMass_.at(species[s]) * thermo_.G(species[s], T);
    }
}


void AFC::MixtureFraction::calculateMeanG
(
    const scalar& H,
    const scalar& S,
    const scalar& T 
)
{
    //- Calculate Gibbs energy with known values
    G_ = thermo_.G(H, S, T);
}
*/


AFC::scalar AFC::MixtureFraction::calculateOmega
(
    const word& species,
    const int Z
) const
{
    //- Concentration map at discrete point
    const map<word, scalar> con = C(Z); 

    return chemistry_.calculateOmega(species, T_[Z], con);
}



/*AFC::scalar AFC::MixtureFraction::calculateHf
(
    const word& species,
    const scalar& T
) const
{
    return thermo_.Hf(species, T);
}*/


void AFC::MixtureFraction::conservation() const
{
    for (int Zi = 0; Zi < nZPoints(); ++Zi)
    {
        conservation(Zi);
    }

    Info<< "\n";
}


void AFC::MixtureFraction::conservation
(
    const int Zi
) const
{
    //- Species
    const wordList& species_ = species();

    scalar Ysum{0};
    scalar Xsum{0};

    forAll(species_, s)
    {
        Ysum += speciesMass_[Zi].at(s); 
        Xsum += speciesMol_[Zi].at(s);
    }

    Info<< "   c-o Mass conservation: " << Zi << " >> " << Ysum << "\n"
        << "       Mole conservation: " << Zi << " >> " << Xsum << "\n";
}


AFC::scalarField AFC::MixtureFraction::Yconservation() const
{
    scalarField Ysum(nZPoints(), 0);

    for (int Zi = 0; Zi < nZPoints(); ++Zi)
    {
        Ysum[Zi] = Yconservation(Zi);
    }

    return Ysum;
}


AFC::scalar AFC::MixtureFraction::Yconservation
(
    const int Zi
) const
{
    //- Species
    const wordList& species_ = species();

    scalar Ysum{0};

    forAll(species_, s)
    {
        Ysum += speciesMass_[Zi].at(s);
    }

    //- Return the absolut value
    return fabs(Ysum);
}

// * * * * * * * * * * * * * * * Update Functions  * * * * * * * * * * * * * //

void AFC::MixtureFraction::updateY
(
    const word& species,
    const scalarField& Y
)
{
    //- Get field access and store new values
    size_t i{0};

    forAll(speciesMass_, Ymap)
    {
        Ymap[species] = Y[i];

        ++i;
    }

    //- Species changed, need to recalculate MMW
    updatedMeanMW(false);
    updatedMeanRho(false);

    actualField_ = "mass";
}


void AFC::MixtureFraction::updateY
(
    const map<word, scalar>& Y_,
    const int Z
)
{
    speciesMass_[Z] = Y_;    
    
    //- Species changed, need to recalculate MMW
    updatedMeanMW(false);
    updatedMeanRho(false);

    actualField_ = "mass";
}


void AFC::MixtureFraction::updateX
(
    const map<word, scalar>& X_,
    const int Z
)
{
    speciesMol_[Z] = X_;    
    
    //- Species changed, need to recalculate MMW
    updatedMeanMW(false);
    updatedMeanRho(false);

    actualField_ = "mol";
}


void AFC::MixtureFraction::updateC
(
    const map<word, scalar>& C_,
    const int Z
)
{
    forAll(C_, c)
    {
        Info<< c.first << " = " << c.second << "\n";
    }
    Info<< "\n";
    forAll(speciesCon_[Z], c)
    {
        Info<< c.first << " = " << c.second << "\n";
    }

    speciesCon_[Z] = C_;    

    //- Species changed, need to recalculate MMW
    updatedMeanMW(false);
    updatedMeanRho(false);

    actualField_ = "con";
}


void AFC::MixtureFraction::updateY()
{
    if (actualField_ == "mass")
    {
        //Info<< "updateY based on mass (return)\n";
        return;
    }

    if (actualField_ == "mol")
    {
        //Info<< "updateY based on mol\n";
        XtoY();
    }
    else if (actualField_ == "con")
    {
        //Info<< "updateY based on con\n";
        CtoY();
    }
}


void AFC::MixtureFraction::updateX()
{
    if (actualField_ == "mol")
    {
        //Info<< "updateX based on mol (return)\n";
        return;
    }

    if (actualField_ == "mass")
    {
        //Info<< "updateX based on mass\n";
        YtoX();
    }
    else if (actualField_ == "con")
    {
        //Info<< "updateX based on con\n";
        CtoX();
    }
}


void AFC::MixtureFraction::updateC()
{
    if (actualField_ == "con")
    {
        //Info<< "updateC based on con (return)\n";
        return;
    }

    if (actualField_ == "mol")
    {
        //Info<< "updateC based on mol\n";
        XtoC();
    }
    else if (actualField_ == "mass")
    {
        //Info<< "updateC based on mass\n";
        YtoC();
    }
}


void AFC::MixtureFraction::updateFields()
{
    updateY();
    updateX();
    updateC();
    updateMeanMW();
    updateRho();
}


void AFC::MixtureFraction::updateT
(
    const scalarField& T
)
{
    T_ = T;
}


void AFC::MixtureFraction::updateRho
(
    const word& method 
)
{
    if (updatedMeanRho())
    {
        return;
    }

    //- Calculate rho using mol fraction X
    if (actualField_ == "mol")
    {
        rhoX();
    }
    else if (actualField_ == "mass")
    {
        rhoY();
    }
    else if (actualField_ == "con")
    {
        rhoC();
    }
    else
    {
        FatalError
        (
            "    Update of mean density can not be done by the method\n"
            "    you choose. Call it with 'mass', 'mol' or 'con' as argument.",
            __FILE__,
            __LINE__
        );
    }

    //- Set to updated
    updatedMeanRho(true);
}


/*void AFC::MixtureFraction::updateCp()
{
    //- Species
    const wordList species = chemistry_.species();

    //- Reset cp
    Cp_ = 0;

    forAll(species, s)
    {
        Cp_ += calculateCp(species[s], temperature_)
            * speciesMol_.at(species[s]); 
    }
}

void AFC::MixtureFraction::updateC()
{
    YtoC();
}*/


void AFC::MixtureFraction::updateOmega()
{
    //- Mixture fraction space (we do not need to calc the boundaries)
    for (int Zi = 1; Zi < nZPoints()-1; ++Zi)
    {
        updateOmega(Zi);
    }
}


void AFC::MixtureFraction::updateOmega
(
    const int Zi 
)
{
    //- Species that are used in combustion
    const wordList& speciesInCombustion = species();

    //- Calculate the source term for each species for all discrete points
    forAll(speciesInCombustion, s)
    {
        updateOmega(Zi, s);
    }
}


void AFC::MixtureFraction::updateOmega
(
    const int Zi,
    const word species
)
{
    //- Calculate the source terms omega for each species and Z
    omega_[Zi][species] = calculateOmega(species, Zi);
}


// * * * * * * * * * * * * * * * * Bool Switch * * * * * * * * * * * * * * * //

bool AFC::MixtureFraction::updatedMeanMW() const
{
    return updatedMeanMW_;
}


void AFC::MixtureFraction::updatedMeanMW
(
    const bool status
)
{
    updatedMeanMW_ = status;
}


bool AFC::MixtureFraction::updatedMeanRho() const
{
    return updatedMeanRho_;
}


void AFC::MixtureFraction::updatedMeanRho
(
    const bool status
)
{
    updatedMeanRho_ = status;
}


void AFC::MixtureFraction::updatedFields
(
    const bool status
)
{
    updatedMeanMW(status);
    updatedMeanRho(status);
}


// * * * * * * * * * * * * Conversation Functions  * * * * * * * * * * * * * //

void AFC::MixtureFraction::YtoX()
{
    //- Update mean molecular weight
    updateMeanMW();

    const wordList& species_ = species();

    //- Mixture fraction space
    for (int Zi = 0; Zi < nZPoints(); ++Zi)
    {

        forAll(species_, s)
        {
            speciesMol_[Zi][s] =
                speciesMass_[Zi].at(s) * MMW_[Zi] / thermo_.MW(s);
        }


        if (debug_)
        {
            scalar sum{0};

            forAll(species_, s)
            {
                sum += speciesMol_[Zi].at(s);
            }

            Info<< "    Z = " << Z(Zi) << ": "
                << "YtoX(), sum of mol = " << sum << endl;
        }
    }
}


void AFC::MixtureFraction::XtoY()
{
    // update mean molecular weight
    updateMeanMW();

    //- Species
    const wordList& species_ = species();

    //- Mixture fraction space
    for (int Zi = 0; Zi < nZPoints(); ++Zi)
    {
        forAll(species_, s)
        {
            speciesMass_[Zi][s] =
                speciesMol_[Zi].at(s) * thermo_.MW(s) / MMW_[Zi];
        }

        if (debug_)
        {
            scalar sum{0};

            forAll(species_, s)
            {
                if (speciesMass_[Zi].at(s) > 0)
                {
                    Info<< s << ": " << speciesMass_[Zi].at(s)
                        << endl;
                }

                sum += speciesMass_[Zi].at(s);
            }

            Info<< "    Z = " << Z(Zi) << ": "
                << "XtoY(), sum of mass = " << sum << endl;
        }
    }
}


void AFC::MixtureFraction::YtoC()
{
    // update mean molecular weight
    updateMeanMW();

    const wordList& species_ = species();

    for (int Zi = 0; Zi < nZPoints(); ++Zi)
    {
        scalar YTbyMW{0};

        forAll(species_, s)
        {
            YTbyMW += speciesMass_[Zi].at(s) * T(Zi) / thermo_.MW(s);
        }

        //- Pressure [Pa]
        const scalar& p = properties_.p();

        //- Unit of concentration [mol/m^3]
        //  Change unit to [mol/cm^3] factor 100cm * 100cm * 100cm / m^3
        forAll(species_, s)
        {
            speciesCon_[Zi][s] =
                (p * speciesMass_[Zi].at(s) / thermo_.MW(s))
                / (AFC::Constants::R * YTbyMW);
        }

        if (debug_)
        {
            scalar sum{0};

            forAll(species_, s)
            {
                sum += speciesCon_[Zi].at(s);
            }

            Info<< "    Z = " << Z(Zi) << ": "
                << "    YtoCon(), sum of concentration = " << sum << endl;
        }
    }
}


void AFC::MixtureFraction::XtoC()
{
    const wordList& species_ = species();

    const scalar& p = properties_.p();

    for (int Zi = 0; Zi < nZPoints(); ++Zi)
    {
        scalar XT{0};

        forAll(species_, s)
        {
            XT += speciesMol_[Zi].at(s) * T(Zi);
        }

        // In [mol/m^3]
        forAll(species_, s)
        {
            speciesCon_[Zi][s] =
                speciesMol_[Zi].at(s) * p / (XT * AFC::Constants::R);
        }
    }
}


void AFC::MixtureFraction::CtoY()
{
    const wordList& species_ = species();

    for (int Zi = 0; Zi < nZPoints(); ++Zi)
    {
        scalar CM{0};

        forAll(species_, s)
        {
            CM += speciesCon_[Zi].at(s) * thermo_.MW(s); 
        }

        forAll(species_, s)
        {
            speciesMass_[Zi][s] = speciesCon_[Zi].at(s) * thermo_.MW(s) / CM;
        }
    }
}


void AFC::MixtureFraction::CtoX()
{
    
}


void AFC::MixtureFraction::rhoX()
{
    //- Update the mean molecular weight
    updateMeanMW();

    const scalar& p = properties_.p();

    for (int Zi = 0; Zi < nZPoints(); ++Zi)
    {
        //- MMW in [g/mol], hence we get [g/m^3], but we save in SI [kg/m^3]
        rho_[Zi] = p * MMW_[Zi] / (AFC::Constants::R * T(Zi)) / scalar(1000);

        if (debug_)
        {
            Info<< "    Mean density (rhoX): " << rho_[Zi] << endl;
        }
    }
}


void AFC::MixtureFraction::rhoY()
{
    //- Update the mean molecular weight
    updateMeanMW();

    const scalar& p = properties_.p();

    for (int Zi = 0; Zi < nZPoints(); ++Zi)
    {
        //- MMW in [g/mol], hence we get [g/m^3], but we save in SI [kg/m^3]
        rho_[Zi] = p * MMW_[Zi] / (AFC::Constants::R * T(Zi)) / scalar(1000);

        if (debug_)
        {
            Info<< "    Mean density (rhoX): " << rho_[Zi] << endl;
        }
    }
}


void AFC::MixtureFraction::rhoC()
{
    const wordList& species = chemistry_.species();

    for (int Zi = 0; Zi < nZPoints(); ++Zi)
    {

        //- Reset rho_
        rho_[Zi] = 0;

        forAll(species, s)
        {
            rho_[Zi] += speciesCon_[Zi].at(s) * thermo_.MW(s);
        }

        if (debug_)
        {
            Info<< "    Mean density (rhoC): " << rho_[Zi] << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

AFC::wordList AFC::MixtureFraction::species() const
{
    return chemistry_.species();
}


AFC::wordList AFC::MixtureFraction::speciesOxidizer() const
{
    return properties_.speciesOxidizer();
}


AFC::wordList AFC::MixtureFraction::speciesFuel() const
{
    return properties_.speciesFuel();
}


int AFC::MixtureFraction::nZPoints() const
{
    return properties_.nZPoints();
}


AFC::scalar AFC::MixtureFraction::Z
(
    const int Z
) const 
{
    return Z_[Z];
}


AFC::scalar AFC::MixtureFraction::T
(
    const int Z 
) const 
{
    return T_[Z];
}


AFC::map<AFC::word, AFC::scalar> AFC::MixtureFraction::X
(
    const int Z 
) const
{
    return speciesMol_[Z];
}


AFC::map<AFC::word, AFC::scalar> AFC::MixtureFraction::Y
(
    const int Z 
) const
{
    return speciesMass_[Z];
}


AFC::map<AFC::word, AFC::scalar> AFC::MixtureFraction::C
(
    const int Z 
) const
{
    return speciesCon_[Z];
}


AFC::map<AFC::word, AFC::scalar> AFC::MixtureFraction::omega
(
    const int Z
) const
{
    return omega_[Z];
}


AFC::scalar AFC::MixtureFraction::sDR() const
{
    return sDR_;
}


AFC::scalarField AFC::MixtureFraction::Z() const 
{
    return Z_;
}


AFC::scalarField AFC::MixtureFraction::rho() const
{
    return rho_;
}


AFC::scalarField AFC::MixtureFraction::T() const
{
    return T_;
}


AFC::scalar AFC::MixtureFraction::p() const
{
    return properties_.p();
}


AFC::List<AFC::map<AFC::word, AFC::scalar> > AFC::MixtureFraction::X() const
{
    return speciesMol_;
}


AFC::List<AFC::map<AFC::word, AFC::scalar> > AFC::MixtureFraction::Y() const
{
    return speciesMass_;
}


AFC::List<AFC::map<AFC::word, AFC::scalar> > AFC::MixtureFraction::C() const
{
    return speciesCon_;
}


AFC::List<AFC::map<AFC::word, AFC::scalar> >
AFC::MixtureFraction::omega() const
{
    return omega_;
}


AFC::scalarField AFC::MixtureFraction::X
(
    const word& species
) const
{
    //- Tmp field
    scalarField XforSpecies(nZPoints(), 0);

    size_t i{0};

    forAll(speciesMol_, discreteX)
    {
        XforSpecies[i] = discreteX.at(species);

        ++i;
    }

    return XforSpecies;
}


AFC::scalarField AFC::MixtureFraction::Y
(
    const word& species
) const
{
    //- Tmp field
    scalarField YforSpecies(nZPoints(), 0);

    size_t i{0};

    forAll(speciesMass_, discreteY)
    {
        YforSpecies[i] = discreteY.at(species);

        ++i;
    }

    return YforSpecies;
}


AFC::scalarField AFC::MixtureFraction::C
(
    const word& species
) const
{
    //- Tmp field
    scalarField CforSpecies(nZPoints(), 0);

    size_t i{0};

    forAll(speciesCon_, discreteC)
    {
        CforSpecies[i] = discreteC.at(species);

        ++i;
    }

    return CforSpecies;
}


AFC::scalarField AFC::MixtureFraction::omega
(
    const word& species 
) const
{
    //- Tmp field
    scalarField omegaforSpecies(nZPoints(), 0);

    size_t i{0};

    forAll(omega_, discreteOmega)
    {
        omegaforSpecies[i] = discreteOmega.at(species);

        ++i;
    }

    return omegaforSpecies;
}


AFC::map<AFC::word, AFC::scalarField> AFC::MixtureFraction::YField() const
{
    //- Tmp map field
    map<word, scalarField> tmp;

    //- Species
    const wordList& speciesI = species();

    //- Loop through all species
    forAll(speciesI, s)
    {
        //- Temporar field that stores Y for species s
        scalarField YforSpecies(nZPoints(), 0);

        size_t Z{0};

        //- Loop through all points
        forAll(speciesMass_, discreteY)
        {
            YforSpecies[Z] = discreteY.at(s);

            ++Z;
        }

        tmp[s] = YforSpecies;
    }

    return tmp;
}


AFC::map<AFC::word, AFC::scalarField> AFC::MixtureFraction::CField() const
{
    //- Tmp map field
    map<word, scalarField> tmp;

    //- Species
    const wordList& speciesI = species();

    //- Loop through all species
    forAll(speciesI, s)
    {
        //- Temporar field that stores [X] for species s
        scalarField CforSpecies(nZPoints(), 0);

        size_t Z{0};

        //- Loop through all points
        forAll(speciesCon_, discreteC)
        {
            CforSpecies[Z] = discreteC.at(s);

            ++Z;
        }

        tmp[s] = CforSpecies;
    }

    return tmp;
}


AFC::map<AFC::word, AFC::scalar> AFC::MixtureFraction::CField
(
    const int Z
) const
{
    //- Tmp map field
    map<word, scalarField> C_ = CField();
    map<word, scalar> CatZ;

    const wordList& speciesI = species();

    forAll(speciesI, s)
    {
        CatZ[s] = C_.at(s)[Z];
    }

    return CatZ;
}


AFC::map<AFC::word, AFC::scalarField> AFC::MixtureFraction::omegaField() const
{
    //- Tmp map field
    map<word, scalarField> tmp;

    //- Species
    const wordList& speciesI = species();

    //- Loop through all species
    forAll(speciesI, s)
    {
        //- Temporar field that stores omega for species s
        scalarField omegaforSpecies(nZPoints(), 0);

        size_t Z{0};

        //- Loop through all points
        forAll(omega_, discreteOmega)
        {
            omegaforSpecies[Z] = discreteOmega.at(s);

            ++Z;
        }

        tmp[s] = omegaforSpecies;
    }

    return tmp;
}


AFC::map<AFC::word, AFC::scalar> AFC::MixtureFraction::omegaField
(
    const int Z
) const
{
    //- Tmp map field
    map<word, scalarField> omega_ = omegaField();
    map<word, scalar> omegaatZ;

    const wordList& speciesI = species();

    forAll(speciesI, s)
    {
        omegaatZ[s] = omega_.at(s)[Z];
    }

    return omegaatZ;
}


AFC::List<int> AFC::MixtureFraction::reacNumbers
(
    const word species
) const
{
    return chemistry_.reacNumbers(species);
}


AFC::wordList AFC::MixtureFraction::speciesInReaction
(
    const int r
) const
{
    return chemistry_.speciesInReaction(r);
}


AFC::wordList AFC::MixtureFraction::speciesProducts
(
    const int r
) const
{
    return chemistry_.speciesProducts(r);
}


AFC::wordList AFC::MixtureFraction::speciesEducts
(
    const int r
) const
{
    return chemistry_.speciesEducts(r);
}


AFC::map<AFC::word, int> AFC::MixtureFraction::nuProducts
(
    const int r
) const
{
    return chemistry_.nuProducts(r);
}


AFC::map<AFC::word, int> AFC::MixtureFraction::nuEducts
(
    const int r
) const
{
    return chemistry_.nuEducts(r);
}


AFC::string AFC::MixtureFraction::elementarReaction
(
    const int r
) const
{
    return chemistry_.elementarReaction(r);
}


AFC::Chemistry AFC::MixtureFraction::chemistry() const
{
    return chemistry_;
}


AFC::Thermo AFC::MixtureFraction::thermo() const
{
    return thermo_;
}


AFC::scalar AFC::MixtureFraction::kf
(
    const int r,
    const int Z
) const
{
    return chemistry_.kf(r, T_[Z]);
}


AFC::scalar AFC::MixtureFraction::kb
(
    const int r,
    const int Z
) const
{
    return chemistry_.kb(r, T_[Z]);
}


// * * * * * * * * * * * * * * * Write output  * * * * * * * * * * * * * * * //

void AFC::MixtureFraction::write() const
{
    Info<< " c-o Save flamelets\n";

    //- Create flamelet folder
    system("mkdir -p flamelet");

    //- Create time folder
    const scalar time = properties_.currentTime();

    //- Command to excecute
    const string cmd = "mkdir -p flamelet/" + toStr(time);

    system(cmd.c_str());

    const string fileName = "flamelet/" + toStr(time) + "/flameletProfil.txt";

    std::filebuf file;

    file.open(fileName, std::ios::out);

    ostream data(&file);

    //- Set scientific notation
    data.setf(std::ios::scientific, std::ios::floatfield);

    //- Build output string and save
    {
        const int nZ = properties_.nZPoints();

        //- Names
        data<< std::setw(6) << "Point"
            << std::setw(15) << "Z"
            << std::setw(15) << "T"
            << std::setw(15) << "rho"
            << std::setw(15) << "mu"
            << std::setw(15) << "lambda"
            << std::setw(15) << "MW" << "  | ";

        //- Species 
        {
            forAll(species(), s)
            {
                data<< std::setw(30) << s << std::setw(15) << "" << "  | ";
            }

            data << "\n";
        }

        //- Fields (mass-fraction, mol-fraction, concentration)
        data<< std::setw(6) << " "
            << std::setw(15) << " "
            << std::setw(15) << " "
            << std::setw(15) << " "
            << std::setw(15) << " "
            << std::setw(15) << " "
            << std::setw(15) << " " << "  | ";

        //- Species
        {
            forAll(species(), s)
            {
                data<< std::setw(15) << "Y"
                    << std::setw(15) << "X"
                    << std::setw(15) << "C" << "  | ";
            }

            data << "\n";
        }

        //- Units
        data<< std::setw(6) << "[-]"
            << std::setw(15) << "[-]"
            << std::setw(15) << "[K]"
            << std::setw(15) << "[kg/m^3]"
            << std::setw(15) << "[kg/m/s]"
            << std::setw(15) << "[W/m/K]"
            << std::setw(15) << "[g/mol]" << "  | ";

        //- Species
        {
            forAll(species(), s)
            {
                data<< std::setw(15) << "[-]"
                    << std::setw(15) << "[-]"
                    << std::setw(15) << "[mol/m^3]" << "  | ";
            }

            data << "\n";
        }

        //- Line
        for (int i = 0; i < 6 * 15 + 6 + 4; ++i)
        {
            data << "=";
        }
        forAll(species(), s)
        {
            for (int i = 0; i < 3 * (15) + 4; ++i)
            {
                data << "=";
            }
        }
        data << "\n";

        for (int i = 0; i < nZ; ++i)
        {
            data<<std::setw(6)<< i+1
                <<std::setw(15)<< Z_[i]
                <<std::setw(15)<< T_[i]
                <<std::setw(15)<< rho_[i]
                <<std::setw(15)<< mu_[i]
                <<std::setw(15)<< lambda_[i]
                <<std::setw(15)<< MMW_[i] << "  | ";

            //- Species mass fraction
            {
                forAll(species(), s)
                {
                    data<< std::setw(15) << speciesMass_[i].at(s)
                        << std::setw(15) << speciesMol_[i].at(s)
                        << std::setw(15) << speciesCon_[i].at(s) << "  | ";
                }

                data << "\n";
            }
        }

        //- Line
        for (int i = 0; i < 6 * 15 + 6 + 4; ++i)
        {
            data << "=";
        }
        forAll(species(), s)
        {
            for (int i = 0; i < 3 * (15) + 4; ++i)
            {
                data << "=";
            }
        }
        data << "\n";
    }

    file.close();
}


AFC::scalar AFC::MixtureFraction::MW
(
    const word species
) const
{
    return thermo_.MW(species);
}


// ************************************************************************* //
