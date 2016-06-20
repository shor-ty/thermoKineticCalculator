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

#include "properties.hpp"
#include "constants.hpp"
#include <cmath>
#include <boost/math/special_functions/erf.hpp>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Properties::Properties
(
    const string& fileName,
    Thermo& thermo,
    const Chemistry& chemistry
)
:
    thermo_(thermo),

    chemistry_(chemistry)
{
    //- Boundarys set all species to zero
    initialBoundary();

    PropertiesReader propertiesReader(fileName);

    //- Init mass and mol fraction map
    {
        const wordList& species = chemistry_.species();

        forAll(species, s)
        {
            fuelX_[species[s]] = 0.;
            fuelY_[species[s]] = 0.;

            oxidizerX_[species[s]] = 0.;
            oxidizerY_[species[s]] = 0.;

            YatZstu_[species[s]] = 0.;
            YatZstb_[species[s]] = 0.;
        }

        const wordList& elements = chemistry_.elements();

        forAll(elements, e)
        {
            oxidizerZj_[elements[e]] = 0.; 
            oxidizerA_[elements[e]] = 0.;

            fuelZj_[elements[e]] = 0.;
            fuelA_[elements[e]] = 0.;

            ZjatZst_[elements[e]] = 0.;
        }
    }

    propertiesReader.read(*this);

    //- Add pressure to thermoData for better handling
    thermo.p(this->p());

    //- Calculate mol or mass fraction of pure streams
    convertFractions();

    //- Check data
    check();

    //- Calculate element mass fraction Zj of fuel and oxidizer
    {
        fuelZj_ = calcZj("F");
        oxidizerZj_ = calcZj("O");
    }

    //- Calculate amount of atoms of fuel and oxidizer
    calcAtomComposition();

    //- Calculate stochiometric mixture fraction
    calcZst();

    //- Calculate strain rate [1/s] at stochiometric condition
    calcStrainRate();

    //- Calculate fuel and oxidizer at stochiometric condition (unburned)
    calcYXatZst();

    //- Calculate element mass fraction Zj at stochiometric condition
    {
        ZjatZst_ = calcZj("ST");
    }

    //- Calc stochiometric coeffs
    calcStochiometricCoeffs(); 

    //- Using element mass fraction we can calculate the species mass fraction
    //  at stochiometric condition only H2O and CO2 + Inerts are available
    //calcYatZstBurned();

    //- Calculate H2O and CO2 mass and mol fraction
    //calcBurnedSpeciesAtZst();

    //- Summary
    //summary();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Properties::~Properties()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AFC::Properties::insertFuel
(
    const word& fuel 
)
{
    fuel_ = fuel;
}


void AFC::Properties::insertOxidizer
(
    const word& oxidizer 
)
{
    oxidizer_ = oxidizer;
}


void AFC::Properties::insertInertGas
(
    const word& inertGas
)
{
    inertGas_ = inertGas;
}


void AFC::Properties::insertMFPoints
(
    const int& mfPoints
)
{
    mfPoints_ = mfPoints;
}


void AFC::Properties::insertVMFPoints
(
    const int& vmfPoints
)
{
    vmfPoints_ = vmfPoints;
}


void AFC::Properties::insertEnthalpyDefects
(
    const scalar& defect
)
{
    defects_.push_back(defect);
}


void AFC::Properties::insertScalarDissipationRates
(
    const scalar& sDR
)
{
    sDRs_.push_back(sDR);
}


void AFC::Properties::insertTemperatureOxidizer
(
    const scalar& TOxidizer
)
{
    TOxidizer_ = TOxidizer; 
}


void AFC::Properties::insertTemperatureFuel
(
    const scalar& TFuel
)
{
    TFuel_ = TFuel;
}


void AFC::Properties::insertCompositionOxidizerMol
(
    const word& species,
    const scalar& molFraction
)
{
    speciesOxidizer_.push_back(species);

    oxidizerX_[species] = molFraction;
}


void AFC::Properties::insertCompositionOxidizerMass
(
    const word& species,
    const scalar& massFraction
)
{
    speciesOxidizer_.push_back(species);

    oxidizerY_[species] = massFraction;
}


void AFC::Properties::insertCompositionFuelMol
(
    const word& species,
    const scalar& molFraction
)
{
    speciesFuel_.push_back(species);

    fuelX_[species] = molFraction;
}


void AFC::Properties::insertCompositionFuelMass
(
    const word& species,
    const scalar& massFraction
)
{
    speciesFuel_.push_back(species);

    fuelY_[species] = massFraction;
}


void AFC::Properties::insertRunTime
(
    const scalar& runTime 
)
{
    runTime_ = runTime;
}


void AFC::Properties::insertWriteControl
(
    const word& writeControl
)
{
    writeControl_ = writeControl;
}


void AFC::Properties::insertWriteControlInterval
(
    const scalar& writeControlInterval
)
{
    writeControlInterval_ = writeControlInterval;
}


void AFC::Properties::insertDeltaT
(
    const scalar& deltaT 
)
{
    deltaT_ = deltaT;
}


void AFC::Properties::insertPressure
(
    const scalar& pressure
)
{
    p_ = pressure;
}


void AFC::Properties::inputMol()
{
    inputMol_ = true;
}


void AFC::Properties::inputMass()
{
    inputMass_ = true;
}


void AFC::Properties::insertInterpreter
(
    const word& keyword
)
{
    interpreter_ = keyword;
}


// * * * * * * * * * * * * * * * Other functions * * * * * * * * * * * * * * //

void AFC::Properties::initialBoundary()
{
    const wordList& species = chemistry_.species();

    forAll(species, s)
    {
        oxidizerX_[species[s]] = 0.;
        oxidizerY_[species[s]] = 0.; 

        fuelX_[species[s]] = 0.;
        fuelY_[species[s]] = 0.;
    }
}


void AFC::Properties::check()
{
    if (mfPoints_ == 0)
    {
        FatalError
        (
            "    No mixtureFractionPoints defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if (vmfPoints_ == 0)
    {
        FatalError
        (
            "    No varianzOfMixtureFractionPoints defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if (TOxidizer_ == 0)
    {
        FatalError
        (
            "    No temperature for oxidizer defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if (TFuel_ == 0)
    {
        FatalError
        (
            "    No temperature for fuel defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if (sDRs_.empty())
    {
        FatalError
        (
            "    No scalar dissipation rates defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if
    (
        oxidizerX_.empty()
     && oxidizerY_.empty()
    )
    {
        FatalError
        (
            "    No oxidizer species defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    if
    (
        fuelX_.empty()
     && fuelY_.empty()
    )
    {
        FatalError
        (
            "    No fuel species defined in the afcDict.",
            __FILE__,
            __LINE__
        );
    }

    scalar sum{0};

    scalar epsilon{1e-5};

    //- Check mass fraction of oxidizer == 1
    {
        //- Check oxidizer stream
        const wordList& species = speciesOxidizer_;

        scalar sumY{0};
        scalar sumX{0};

        forAll(species, s)
        {
            sumY += oxidizerY_.at(species[s]);
            sumX += oxidizerX_.at(species[s]);
        }

        if (std::abs(1-sumY) > epsilon || std::abs(1-sumX) > epsilon)
        {
            const string valueY = std::to_string(sumY);
            const string valueX = std::to_string(sumX);

            FatalError
            (
                "    Mass/mol fraction of oxidizer is not 1.\n" 
                "    Value of Y: " + valueY + " and X: " + valueX,
                __FILE__,
                __LINE__
            );
        }
    }

    //- Check mass and mol fraction of fuel == 1
    {
        //- Check oxidizer stream
        const wordList& species = speciesFuel_;

        scalar sumY{0};
        scalar sumX{0};

        forAll(species, s)
        {
            sumY += fuelY_.at(species[s]);
            sumX += fuelX_.at(species[s]);
        }

        if (std::abs(1-sumY) > epsilon || std::abs(1-sumX) > epsilon)
        {
            const string valueY = std::to_string(sumY);
            const string valueX = std::to_string(sumX);

            FatalError
            (
                "    Mass/mol fraction of fuel is not 1.\n" 
                "    Value of Y: " + valueY + " and X: " + valueX,
                __FILE__,
                __LINE__
            );
        }
    }

    //- Check if pressure is set
    if (p_ <= 0)
    {
        FatalError
        (
            "    Pressure is not set or set not correct.\n"
            "    Please check the pressure in afcDict.",
            __FILE__,
            __LINE__
        );
    }

    //- Check if fuel is set
    if (fuel().empty())
    {
        FatalError
        (
            "    No fuel species is set.\n"
            "    Please check the fuel keyword in afcDict.",
            __FILE__,
            __LINE__
        );
    }

    //- Check if oxidizer is set
    if (oxidizer().empty())
    {
        FatalError
        (
            "    No oxidizer species is set.\n"
            "    Please check the oxidizer keyword in afcDict.",
            __FILE__,
            __LINE__
        );
    }
}


void AFC::Properties::convertFractions()
{
    if (inputMol_) 
    {
        //- Convert to mass fraction
        XtoY();
    }
    else if (inputMass_)
    {
        //- Convert to mol fraction
        YtoX();
    }
}


void AFC::Properties::XtoY()
{

    //- Oxidizer stream
    {
        const wordList& species = speciesOxidizer_;

        //- Mean molecular weight
        scalar MMW{0};

        forAll(species, s)
        {
            MMW += oxidizerX_.at(species[s]) * thermo_.MW(species[s]);
        }

        forAll(species, s)
        {
            oxidizerY_[species[s]] =
                oxidizerX_.at(species[s]) * thermo_.MW(species[s]) / MMW;
        }
    }

    //- Fuel stream
    {
        const wordList& species = speciesFuel_;

        //- Mean molecular weight
        scalar MMW{0};

        forAll(species, s)
        {
            MMW += fuelX_.at(species[s]) * thermo_.MW(species[s]);
        }

        forAll(species, s)
        {
            fuelY_[species[s]] =
                fuelX_.at(species[s]) * thermo_.MW(species[s]) / MMW;
        }
    }
}


void AFC::Properties::YtoX()
{
    //- Oxidizer stream
    {
        const wordList& species = speciesOxidizer_;

        //- Mean molecular weight
        scalar denominator{0};

        forAll(species, s)
        {
            denominator += oxidizerY_.at(species[s]) / thermo_.MW(species[s]);
        }

        const scalar& MMW = 1 / denominator;

        forAll(species, s)
        {
            oxidizerX_[species[s]] =
                oxidizerY_.at(species[s]) * MMW / thermo_.MW(species[s]); 
        }
    }
    
    //- Fuel stream
    {
        const wordList& species = speciesFuel_;

        //- Mean molecular weight
        scalar denominator{0};

        forAll(species, s)
        {
            denominator += fuelY_.at(species[s]) / thermo_.MW(species[s]);
        }

        const scalar& MMW = 1 / denominator;

        forAll(species, s)
        {
            fuelX_[species[s]] =
                fuelY_.at(species[s]) * MMW / thermo_.MW(species[s]);
        }
    }
}


void AFC::Properties::HFuelAdiabatic()
{
    //- Fuel species
    const wordList& species = speciesFuel();

    //- Fuel mass composition
    const map<word, scalar>& massFraction = fuelY();
    const map<word, scalar>& molFraction = fuelX();

    //- Temperature of fuel
    const scalar& T = fuelTemperature();

    //- Adiabatic enthalpy
    scalar tmp{0};

    forAll(species, s)
    {
        //- Enthalpy for species s [J/mol]
        const scalar& H = thermo_.H(species[s], T);

        //- Molecular weight [g/mol]
        const scalar& MW = thermo_.MW(species[s]);

        //- H_mol / MW * X
        tmp += H / MW * massFraction.at(species[s]) * 1000;
    }

    fuelH_ = tmp;
}


void AFC::Properties::HOxidizerAdiabatic()
{
    //- Oxidizer species
    const wordList& species = speciesOxidizer();

    //- Oxidizer mass composition
    const map<word, scalar>& massFraction = oxidizerY();

    //- Temperature of fuel
    const scalar& T = oxidizerTemperature();

    //- Adiabatic enthalpy
    scalar tmp{0};

    forAll(species, s)
    {
        //- Enthalpy for species s [J/mol]
        const scalar& H = thermo_.H(species[s], T);

        //- Molecular weight [g/mol]
        const scalar& MW = thermo_.MW(species[s]);

        tmp += H / MW * massFraction.at(species[s]) * 1000;
    }

    oxidizerH_ = tmp;
}


void AFC::Properties::calcProperties()
{
    //- Calc adiabatic fuel enthalpy
    HFuelAdiabatic(); 

    //- Calc adibatic oxidizer enthalpy
    HOxidizerAdiabatic();

}


AFC::map<AFC::word, AFC::scalar> AFC::Properties::calcZj
(
    const word& OF 
)
{
    const wordList& species = chemistry_.species();

    //- Tmp Zj
    map<word, scalar> ZjMap;

    //- Elements used
    const wordList& elements = chemistry_.elements();

    //- Element loop
    forAll(elements, e)
    {
        //- Element molecular weight (atom weight)
        const scalar& MWElement = AFC::Constants::AW.at(elements[e]);

        //- Temp Zj 
        scalar Zj{0};

        //- Loop through species
        forAll(species, s)
        {
            //- All elements in species
            const wordList& elementsInSpecies =
                thermo_.elementsInSpecies(species[s]);

            //- Number of elements in species
            const scalarList& nElem = thermo_.elementsFactors(species[s]);
    
            //- Find index of element in species and get atom number
            scalar a{0};

            forAll(elementsInSpecies, s)
            {
                if (elementsInSpecies[s] == elements[e])
                {
                    a = nElem[s];
                }
            }

            const scalar& MWSpecies = thermo_.MW(species[s]);

            scalar Y{0};

            if(OF == "O")
            {
                Y = oxidizerY_.at(species[s]);
            }
            else if (OF == "F")
            {
                Y = fuelY_.at(species[s]);
            }
            else if (OF == "ST")
            {
                Y = YatZstu_.at(species[s]);
            }

            Zj += a * MWElement / MWSpecies * Y ;
        }

        ZjMap[elements[e]] = Zj;
    }

    return ZjMap;
}


void AFC::Properties::calcZst()
{
    //- Moleculare weight of atoms
    const scalar& MWC = AFC::Constants::AW.at("C");
    const scalar& MWH = AFC::Constants::AW.at("H");
    const scalar& MWO = AFC::Constants::AW.at("O");

    //- Ratios (to burn 1 C or 1 H)
    const scalar ominC = 2 * MWO / MWC;
    const scalar ominH = 0.5 * MWO / MWH;

    //- Element mass fraction of fuel
    const map<word, scalar>& Zj_F = fuelZj();
    const scalar& Zj_H = Zj_F.at("H");
    const scalar& Zj_C = Zj_F.at("C");
    const scalar& Zj_O = Zj_F.at("O");

    //- omin == nu for stochiometric
    const scalar omin = ominC * Zj_C + ominH * Zj_H - Zj_O;

    omin_ = omin;

    //- Mass fraction of O2 in oxidizer
    const scalar& YO_O2 = oxidizerY("O2");

    //- Mass fraction of fuel in fuel
    const wordList& fuelSpecies = speciesFuel();

    scalar YF_F{0};

    forAll(fuelSpecies, s)
    {
        //- If C or H is in species == fuel
        //  TODO - what if  CO2 or H2O = ? 
        const wordList& atoms = thermo_.elementsInSpecies(fuelSpecies[s]);

        bool fuel{false};

        forAll(atoms, a)
        {
            if(atoms[a] == "C" || atoms[a] == "H")
            {
                fuel = true;
            }
        }

        if (fuel)
        {
            YF_F += fuelY(fuelSpecies[s]);
        }
    }

    //- Save 
    Zst_ = pow(1 + (omin * YF_F / YO_O2), -1);
}


void AFC::Properties::calcStrainRate()
{
    //- Calculate strain rate with stochiometric mixture fraction Zst
    //  and stochiometric scalar dissipation rate
    {
        //- Scalar dissipation rates
        const scalarField& chis = sDRs();
        const scalar Z_st = Zst();

        forAll(chis, c)
        {
            as_.push_back
            (
                chis[c] * M_PI / exp(-2*pow(boost::math::erfc_inv(2*Z_st), 2)) 
            );
        }
    }
}


void AFC::Properties::calcYXatZst()
{
    //- Species of fuel and oxidizer 
    const wordList& speciesF = speciesFuel();
    const wordList& speciesO = speciesOxidizer();

    //- Oxidizer mass fraction
    const map<word, scalar>& oxidizerMassFraction = oxidizerY();
    
    //- Fuel mass fraction
    const map<word, scalar>& fuelMassFraction = fuelY();

    //- Stochiometric mixture fraction Zst
    const scalar& Z_st = Zst();

    //- All species 
    const wordList& species = chemistry_.species();

    //- Mass fraction
    forAll(species, i)
    {
        if
        (
            oxidizerMassFraction.at(species[i]) > 1e-10
         || fuelMassFraction.at(species[i]) > 1e-10
        )
        {
            YatZstu_[species[i]] =
                oxidizerMassFraction.at(species[i]) * (1 - Z_st)
              + fuelMassFraction.at(species[i]) * Z_st;
        }
        else
        {
            YatZstu_[species[i]] = 0;
        }
    }
    //- Mol fraction 
    {
        //- Mean molecular weight
        scalar denominator{0};

        forAll(species, s)
        {
            denominator += YatZstu_.at(species[s]) / thermo_.MW(species[s]);
        }

        const scalar& MMW = 1 / denominator;

        forAll(species, s)
        {
            XatZstu_[species[s]] =
                YatZstu_.at(species[s]) * MMW / thermo_.MW(species[s]); 
        }
    }
}


void AFC::Properties::calcAdiabaticTemperature()
{
    //- Very simple first approach for the flame temperature at Zst
    //  Combustion is calculated using CnHm + O2 = H2O + CO2
    
    //- Fuel and oxidizer species
    const word& F = fuel();
    const word& O = oxidizer();

    //- Get amount of H and C of fuel
    const wordList& elements = thermo_.elementsInSpecies(F);
    const scalarList& nElements = thermo_.elementsFactors(F);

    //- Get amount of oxigen, H2O and CO2
    scalar nOxigen{0};
    scalar nH2O{0};
    scalar nCO2{0};

    forAll(elements, e)
    {
        if (elements[e] == "H")
        {
            //- For H we need 1/2 O (H2O)
            nOxigen += nElements[e] / 2.;
            nH2O += nElements[e] / 2.;
        }
        else if (elements[e] == "C")
        {
            //- For C we need 2 O (CO2)
            nOxigen += nElements[e] * 2.;
            nCO2 += nElements[e] * 1.;
        }
        else if (elements[e] == "O")
        {
            //- If we have O in fuel 
            nOxigen -= nElements[e];
        }
    }

    //- The amount of O2
    nOxigen /= 2.;

    //- For burning 1g of F we need nOxigen g, therefore the mass
    //  at unburned state is the sum:
    const scalar m_u = 1. + nOxigen;
    const scalar YF = 1. / m_u;
    const scalar YO = nOxigen / m_u;

    //- Hence we end up with nCO2 and nH2O
    const scalar m_b = nH2O + nCO2;
    const scalar YCO2 = nCO2 / m_b;
    const scalar YH2O = nH2O / m_b;

    //- Unburned and two temperature at Zst [K]
    const scalar& TFuel = fuelTemperature();
    const scalar& TOxid = oxidizerTemperature();
    const scalar Tu = (TFuel - TOxid) * Zst() + TFuel;

    //- T1 high temperature, T2 low temperature, Tm mean of T1 and T2
    //  T1 > Tb > T2
    scalar T1{5000};
    scalar T2{500};

    //- Molar enthalpy converted to specific of unburned state
    const scalar& hF = thermo_.H(F, Tu) / thermo_.MW(F);
    const scalar& hO2 = thermo_.H(O, Tu) / thermo_.MW(O);
    const scalar hu = YF * hF + YO * hO2;

    unsigned int iter{500};

    for(unsigned int i=0; i<iter; i++)
    {
        scalar Tm = (T1 + T2)/2;
        
        //- Molar enthalpy converted to specific for burned state T1, T2, Tm
        const scalar& hCO2_T1 = thermo_.H("CO2", T1) / thermo_.MW("CO2");
        const scalar& hH2O_T1 = thermo_.H("H2O", T1) / thermo_.MW("H2O");

        const scalar& hCO2_T2 = thermo_.H("CO2", T2) / thermo_.MW("CO2");
        const scalar& hH2O_T2 = thermo_.H("H2O", T2) / thermo_.MW("H2O");

        const scalar& hCO2_Tm = thermo_.H("CO2", Tm) / thermo_.MW("CO2");
        const scalar& hH2O_Tm = thermo_.H("H2O", Tm) / thermo_.MW("H2O");

        //- Specific enthalpy of burned state
        const scalar hb_T1 = YCO2 * hCO2_T1 + YH2O * hH2O_T1;
        const scalar hb_T2 = YCO2 * hCO2_T2 + YH2O * hH2O_T2;
        const scalar hb_Tm = YCO2 * hCO2_Tm + YH2O * hH2O_Tm;

        //- Mean molar heat capacity converted to specific for burned state
        const scalar T1m = (T1 + Tu)/2.;
        const scalar T2m = (T2 + Tu)/2.;
        const scalar Tmm = (Tm + Tu)/2.;
        
        const scalar& cpCO2_T1 = thermo_.cp("CO2", T1m) / thermo_.MW("CO2");
        const scalar& cpH2O_T1 = thermo_.cp("H2O", T1m) / thermo_.MW("H2O");
        
        const scalar& cpCO2_T2 = thermo_.cp("CO2", T2m) / thermo_.MW("CO2");
        const scalar& cpH2O_T2 = thermo_.cp("H2O", T2m) / thermo_.MW("H2O");

        const scalar& cpCO2_Tm = thermo_.cp("CO2", Tmm) / thermo_.MW("CO2");
        const scalar& cpH2O_Tm = thermo_.cp("H2O", Tmm) / thermo_.MW("H2O");

        //- Mean specific heat capacity of burned state
        const scalar cp_T1 = YCO2 * cpCO2_T1 + YH2O * cpH2O_T1;
        const scalar cp_T2 = YCO2 * cpCO2_T2 + YH2O * cpH2O_T2;
        const scalar cp_Tm = YCO2 * cpCO2_Tm + YH2O * cpH2O_Tm;

        //- Adiabatic system, without dissociation, simplified without inert gas
        //  and other species. Could be improved, but is only one guess
        //  hu + dQ = hb; dQ = cpb * (Tb - Tu) 
        //  Better calculation would be that one used by Peters (script) 
        //  Here, hu has to be equal to hb - dQ
        const scalar dQ1 = cp_T1 * (T1 - Tu);
        const scalar dQ2 = cp_T2 * (T2 - Tu);
        const scalar dQm = cp_Tm * (Tm - Tu);

        const scalar hb1 = hb_T1 + dQ1;
        const scalar hb2 = hb_T2 + dQ2;
        const scalar hbm = hb_Tm + dQm;

        if (hb1 > hu && hu > hbm)
        {
            //- Change T2 to Tm
            T2 = Tm;
        }
        else if (hbm > hu && hu > hb2)
        {
            //- Change T1 to Tm
            T1 = Tm;
        }
        else
        {
            //- Some problem in calculation
            FatalError
            (
                "    Problems during estimating the adiabatic flame"
                " temperature.\n    The enthalpy comparison has problems.",
                __FILE__,
                __LINE__
            );
        }

        if(abs(hu - hb1) < 1e-5)
        {
            //- Some assumption nothing that is valid
            /*if (Tm > 2800 && Tm < 3200)
            {
                Tm *= 0.75;
            }
            else if (Tm >= 3200)
            {
                Tm *= 0.75;
            }*/

            Tadiabatic_ = Tm;

            break;
        }
        
        //- if too many iterations error
        if (i == iter)
        {
            FatalError
            (
                "    The adiabatic flame temperature exceed the maximum\n"
                "    iteration of " + std::to_string(iter) + ".",
                __FILE__,
                __LINE__
            );
        }
    }
}


void AFC::Properties::calcAtomComposition()
{

    //- Oxidizer
    {
        //- a) Get all species of oxidizer
        const wordList& species = speciesOxidizer();

        //- b) Get all elements that are in oxidizer
        wordList elements;
        scalarList n;

        forAll(species, s)
        {
            const wordList& atoms = thermo_.elementsInSpecies(species[s]);
            const scalarList& nAtoms = thermo_.elementsFactors(species[s]);

            forAll(atoms, a)
            {
                elements.push_back(atoms[a]);
                n.push_back(nAtoms[a]);
            }
        }

        //- c) get global amount of atoms in oxidizer
        wordList tmp;
        scalarList nTmp;

        forAll(elements, e)
        {
            bool found{false};

            forAll(tmp, t)
            {
                if(tmp[t] == elements[e])
                {
                    //- Add number to atom
                    nTmp[t] += n[e];
                    found = true;
                    break;
                }
            }

            //- First loop or if atom not yet insert
            if (tmp.empty() || ! found)
            {
                tmp.push_back(elements[e]);
                nTmp.push_back(n[e]);
            }
        }

        //- d) set atom composition 
        forAll(tmp, t)
        {
            oxidizerA_[tmp[t]] = nTmp[t];
        }
    }

    //- Fuel
    {
        //- a) Get all species of fuel
        const wordList& species = speciesFuel();

        //- b) Get all elements that are in fuel 
        wordList elements;
        scalarList n;

        forAll(species, s)
        {
            const wordList& atoms = thermo_.elementsInSpecies(species[s]);
            const scalarList& nAtoms = thermo_.elementsFactors(species[s]);

            forAll(atoms, a)
            {
                elements.push_back(atoms[a]);
                n.push_back(nAtoms[a]);
            }
        }

        //- c) get global amount of atoms in fuel
        wordList tmp;
        scalarList nTmp;

        forAll(elements, e)
        {
            bool found{false};

            forAll(tmp, t)
            {
                if(tmp[t] == elements[e])
                {
                    //- Add number to atom
                    nTmp[t] += n[e];
                    found = true;
                    break;
                }
            }

            //- First loop or if atom not yet insert
            if (tmp.empty() || ! found)
            {
                tmp.push_back(elements[e]);
                nTmp.push_back(n[e]);
            }
        }

        //- d) set atom composition
        forAll(tmp, t)
        {
            fuelA_[tmp[t]] = nTmp[t];
        }
    }
}


void AFC::Properties::calcStochiometricCoeffs()
{
    //- Assumption (Warnatz)
    const int nuB = 1;

    nuCO2_ = nuB * fuelA_.at("C");
    nuH2O_ = nuB * fuelA_.at("H") / 2;
    nuO2_ = nuCO2_ + nuH2O_ / 2 - nuB * fuelA_.at("O") / 2;
}


void AFC::Properties::calcBurnedSpeciesAtZst()
{
    //- Mol of products
    scalar sumMol = nuCO2_ + nuH2O_;

    //- Mol fraction of CO2 and H2O 
    XatZstb_["H2O"] = nuH2O_ / sumMol;
    XatZstb_["CO2"] = nuCO2_ / sumMol;

    const wordList& species{"H2O", "CO2"};

    //- Mean molecular weight
    scalar MMW{0};

    MMW += XatZstb_.at("H2O") * thermo_.MW("H2O");
    MMW += XatZstb_.at("CO2") * thermo_.MW("CO2");

    //- Calculate mass fraction 
    YatZstb_["H2O"] = XatZstb_.at("H2O") * thermo_.MW("H2O") / MMW;
    YatZstb_["CO2"] = XatZstb_.at("CO2") * thermo_.MW("CO2") / MMW;
}


// * * * * * * * * * * * * * * * Return Functions  * * * * * * * * * * * * * //

AFC::word AFC::Properties::fuel() const
{
    return fuel_;
}


AFC::word AFC::Properties::oxidizer() const
{
    return oxidizer_;
}


AFC::word AFC::Properties::inertGas() const
{
    return inertGas_;
}


AFC::wordList AFC::Properties::speciesOxidizer() const
{
    return speciesOxidizer_;
}


AFC::map<AFC::word, AFC::scalar> AFC::Properties::oxidizerX() const
{
    return oxidizerX_;
}


AFC::map<AFC::word, AFC::scalar> AFC::Properties::oxidizerY() const
{
    return oxidizerY_;
}


AFC::scalar AFC::Properties::oxidizerY
(
    const word& species 
) const
{
    return oxidizerY_.at(species);
}


AFC::map<AFC::word, AFC::scalar> AFC::Properties::oxidizerZj() const
{
    return oxidizerZj_;
}


AFC::wordList AFC::Properties::speciesFuel() const
{
    return speciesFuel_;
}


AFC::map<AFC::word, AFC::scalar> AFC::Properties::fuelX() const
{
    return fuelX_;
}


AFC::map<AFC::word, AFC::scalar> AFC::Properties::fuelY() const
{
    return fuelY_;
}


AFC::scalar AFC::Properties::fuelY
(
    const word& species
) const
{
    return fuelY_.at(species);
}


AFC::map<AFC::word, AFC::scalar> AFC::Properties::fuelZj() const
{
    return fuelZj_;
}


AFC::scalarField AFC::Properties::defects() const
{
    return defects_;
}


AFC::scalarField AFC::Properties::sDRs() const
{
    return sDRs_;
}


AFC::scalar AFC::Properties::sDRs
(
    const int& i
) const
{
    return sDRs_[i];
}


int AFC::Properties::mfPoints() const
{
    return mfPoints_;
}


int AFC::Properties::vmfPoints() const
{
    return vmfPoints_;
}


AFC::scalar AFC::Properties::oxidizerTemperature() const
{
    return TOxidizer_;
}


AFC::scalar AFC::Properties::fuelTemperature() const
{
    return TFuel_;
}


AFC::scalar AFC::Properties::runTime() const
{
    return runTime_;
}


AFC::scalar AFC::Properties::deltaT() const
{
    return deltaT_;
}


unsigned int AFC::Properties::nDefects() const
{
    return defects_.size();
}


AFC::scalar AFC::Properties::defect
(
    const int& defectNo
) const
{
    return defects_[defectNo];
}


AFC::scalar AFC::Properties::p() const
{
    return p_;
}


AFC::word AFC::Properties::input() const
{
    if
    (
        !inputMol_
     && !inputMass_
    )
    {
        FatalError
        (
            "    Input problems for mass or mol fraction",
            __FILE__,
            __LINE__
        );
    }

    word ret{"none"};

    if (inputMol_)
    {
        ret = "mol";
    }
    else if (inputMass_)
    {
        ret = "mass";
    }

    return ret;
}


AFC::word AFC::Properties::interpreter() const
{
    return interpreter_;
}


AFC::scalar AFC::Properties::Zst() const
{
    return Zst_;
}


AFC::scalar AFC::Properties::YatZstu
(
    const word& species
) const
{
    return YatZstu_.at(species);
}


AFC::scalar AFC::Properties::YatZstb
(
    const word& species
) const
{
    return YatZstb_.at(species);
}


AFC::scalar AFC::Properties::Tadiabatic() const
{
    return Tadiabatic_;
}


// ************************************************************************* //
