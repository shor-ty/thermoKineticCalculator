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

#include "typedef.hpp"
#include "jacobian.hpp"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AFC::Jacobian::Jacobian
(
    const Chemistry& chem 
)
:
    chem_(chem)
{
    if (debug_)
    {
        Info<< "Constructor Jacobian \n" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AFC::Jacobian::~Jacobian()
{
    if (debug_)
    {
        Info<< "Destructor Jacobian \n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

void AFC::Jacobian::jacobian
(
    const scalar T,
    const scalar p,
    const scalar t,
    const map<word, scalar>& con,
    map<word, scalar>& dcdt,
    Matrix& dcdc
)
{
    if (debug_)
    {
        Info<< " --> Enter AFC::Jacobian::jacobian" << endl;
    }

    //- Calculate the Jacobian matrix

    //- Species
    const wordList& species = chem_.species();

    //- Number of species
    const size_t nSpecies = species.size();

    //- Copy of concentration
    map<word, scalar> con1 = con;

    //- Make sure that c1 is positive)
    forAll(species, s)
    {
        con1[s] = max(con1.at(s), scalar(0));
    }

    //- Reset the Jacobian matrix (set values to zero)
    dcdc.reset();

    //- Calculate the rate of each species
    dcdt = chem_.calculateOmega(T, con1);

    //- Build Jacobian matrix (n x n)
    Matrix& J = dcdc;

    //  The first column is the derivation with repect to the first species
    //  i = row
    //  j = column
    
    //- Row 
    for (size_t i = 0; i<nSpecies; ++i)
    {
        //- f([x]) :: X -> Species1
        const word species1 = species[i];

        //- Get information about the reactions where species1 is involved
        const List<int>& inReaction = chem_.reacNumbers(species1);

        //- Column 
        for (size_t j = 0; j<nSpecies; ++j)
        {
            //- The species2 is the derivation guy
            const word species2 = species[j];

            Info<< "Jacobian for " << species1 << " | " << species2 << endl;

            //- Entry in the Jacobian matrix at position i j
            scalar Jvalue{0};

            //- Loop through all reactions
            for (size_t ri = 0; ri<inReaction.size(); ++ri)
            {
                //- Reaction number
                const size_t r = inReaction[ri];

                //- Get the species that are in the reaction r
                const wordList& speciesInReac = chem_.speciesInReaction(r);

                //- Check if species2 is in the reaction
                bool cont{false};

                forAll(speciesInReac, s)
                {
                    if (s == species2)
                    {
                        cont = true;
                        break;
                    }
                }

                //- If species is in the reaction, we need to make the
                //  derivation. If it is not included, the stuff gets zero
                //  and we do not have to take care about that
                if (cont)
                {
                    //- Get educt/product species and the corresponding 
                    //  stochiometetric factors
                    const List<word>& eductSpecies = chem_.speciesEducts(r);
                    const List<word>& productSpecies = chem_.speciesProducts(r);
                    const map<word, int> nuEduct = chem_.nuEducts(r);
                    const map<word, int> nuProduct = chem_.nuProducts(r);

                    //- TODO other workaround to expensive
                    const scalar kf = chem_.kf(r, T);
                    const scalar kb = chem_.kb(r, T);

                    Info<< "  Investigating into " << chem_.elementarReaction(r) << "\n";

                    //- Return the value of the derivation of the elementar
                    //  reaction with respect to species2
                    Jvalue += derivationOfReaction
                        (
                            species1,
                            species2,
                            eductSpecies,
                            productSpecies,
                            nuEduct,
                            nuProduct,
                            kf,
                            kb,
                            con
                        );
                }

            //- Loop reactions
            }

            //- Jvalue calculated - now set it
            J(i,j, Jvalue);

        //- Loop derivative species
        }

        Info<< "\n======================================\n";

    //- Loop function
    }

    if (debug_)
    {
        Info<< " --> Leave AFC::Jacobian::jacobian" << endl;
    }
}


AFC::scalar AFC::Jacobian::derivationOfReaction
(
    const word species,
    const word speciesOfDerivation,
    const wordList& eductSpecies,
    const wordList& productSpecies,
    const map<word, int>& nuEduct,
    const map<word, int>& nuProduct,
    const scalar kf,
    const scalar kb,
    const map<word, scalar>& con
) const
{
    if (debug_)
    {
        Info<< " --> Enter AFC::Jacobian::derivationOfReaction" << endl;
    }

    //- This function calculates the derivation of the elementar reaction
    //  with respect to the speciesOfDerivation
   
    //- If the concentration of speciesOfDerivation is zero, -> return 0
    if (con.at(speciesOfDerivation) < 1e-15)
    {
        return 0;
    }

    
    //- Check where the species is included
    bool inEduct{false};
    bool inProduct{false};

    forAll(eductSpecies, s)
    {
        if (s == speciesOfDerivation)
        {
            inEduct = true;
            break;
        }
    }

    forAll(productSpecies, s)
    {
        if (s == speciesOfDerivation)
        {
            inProduct = true;
            break;
        }
    }

    //- Calculate the product of the concentrations
    //  TODO [M] not included
    scalar conEduct{1};
    scalar conProduct{1};
        
    if (inEduct)
    {
        //- Linear 
        if (abs(nuEduct.at(speciesOfDerivation)) == 1)
        {
            forAll(eductSpecies, s)
            {
                if (s != speciesOfDerivation)
                {
                    conEduct *= pow(con.at(s), abs(nuEduct.at(s)));
                }
            }
        }
        //- Non-Linear (Power x)
        else if (abs(nuEduct.at(speciesOfDerivation)) > 1)
        {
            forAll(eductSpecies, s)
            {
                if (s != speciesOfDerivation)
                {
                    conEduct *= pow(con.at(s), abs(nuEduct.at(s)));
                }
                else
                {
                    conEduct *= abs(nuEduct.at(s))
                        * pow(con.at(s), abs(nuEduct.at(s))-1);
                }
            }
        }
        else
        {
            FatalError
            (
                "    The stochiometric factor of species "
                + speciesOfDerivation + " is not 1 or 2.\n"
                " Not implemented | or reactions wrong.",
                __FILE__,
                __LINE__
            );
        }
    }
    else
    {
        conEduct = 0;
    }


    if (inProduct)
    {
        
        //- Linear 
        if (nuProduct.at(speciesOfDerivation) == 1)
        {
            forAll(productSpecies, s)
            {
                if (s != speciesOfDerivation)
                {
                    conProduct *= pow(con.at(s), abs(nuProduct.at(s)));
                }
            }
        }
        //- Non-Linear (Power x)
        else if (nuProduct.at(speciesOfDerivation) > 1)
        {
            forAll(productSpecies, s)
            {
                if (s != speciesOfDerivation)
                {
                    conProduct *= pow(con.at(s), abs(nuProduct.at(s)));
                }
                else
                {
                    conProduct *= abs(nuProduct.at(s))
                        * pow(con.at(s), abs(nuProduct.at(s))-1);
                }
            }
        }
        else
        {
            FatalError
            (
                "    The stochiometric factor of species "
                + speciesOfDerivation + " is not 1 or 2.\n"
                " Not implemented | or reactions wrong.",
                __FILE__,
                __LINE__
            );
        }
    }
    else
    {
        conProduct = 0;
    }

    /*if (inEduct)
    {
        Info<< "    --- in educt site\n";
    }
    if (inProduct)
    {
        Info<< "    --- in product site\n";
    }*/

    //- Get pre-factor (nu'' - nu') :: based on the fact that we already
    //  know the right value, we just have to check if the species
    //  is within the product or educt side. We need + because educt 
    //  took already a minus sign.
    scalar nuSpecies{0};

    if (nuEduct.count(species) && !nuProduct.count(species))
    {
        nuSpecies = nuEduct.at(species);
    }
    else if (!nuEduct.count(species) && nuProduct.count(species))
    {
        nuSpecies = nuProduct.at(species);
    }
    else if (nuEduct.count(species) && nuProduct.count(species))
    {
        nuSpecies = nuProduct.at(species) + nuEduct.at(species);
    }
    else
    {
        FatalError
        (
            "Not implemented. Error.",
            __FILE__,
            __LINE__
        );
    }

    /*
    {

    Info<< std::setw(20) << "nu'' - nu'  = " <<std::setw(20) <<  nuSpecies << "\n";
    Info<< std::setw(20) << "kf = " <<std::setw(20) <<  kf << "\n";
    Info<< std::setw(20) << "kb = " <<std::setw(20) <<  kb << "\n";
    Info<< std::setw(20) << " conEduc = "<< std::setw(20) << conEduct << "\n";
    Info<< std::setw(20) << " conProd = "<< std::setw(20) << conProduct<< "\n";

    //- Return the derivation value
    Info << std::setw(20)<< " === " << std::setw(20)
        <<  nuSpecies * (kf * conEduct - kb * conProduct) << "\n\n";
    }
    */

    if (debug_)
    {
        Info<< " --> Leave AFC::Jacobian::derivationOfReaction" << endl;
    }

    return nuSpecies * (kf * conEduct - kb * conProduct);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
