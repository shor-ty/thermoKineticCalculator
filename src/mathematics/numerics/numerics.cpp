/*---------------------------------------------------------------------------*\
  c-o-o-c-o-o-o             |
  |     |     A utomatic    | Open Source Flamelet
  c-o-o-c     F lamelet     | 
  |     |     C onstructor  | Copyright (C) 2020 Holzmann CFD
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

#include "definitions.hpp"
#include "numerics.hpp"
//#include "ODE.hpp"
#include "rosenbrock.hpp"
//#include "euler.hpp"
#include <math.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
AFC::Numerics<Type>::Numerics
(
    Chemistry& chem 
)
{
    if (debug_)
    {
        Info<< "Constructor Numerics \n" << endl;
    }

    //- Create new object
    //    std::unique_ptr<Type> foo= std::make_unique<Type>(chem);
    //ode_ = std::unique_ptr<Type>();
    ode_ = new Type(chem);

    //- TODO readable
    deltaTChem_ = scalarField(11, 1e-7);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
AFC::Numerics<Type>::~Numerics()
{
    if (debug_)
    {
        Info<< "Destructor Numerics \n" << endl;
    }

    //- Delete the ode_ pointer to free the memory
    delete ode_;

}


// * * * * * * * * * * * * * * * Member function * * * * * * * * * * * * * * //

template<class Type>
void AFC::Numerics<Type>::solveForInitialSolution
(
    MixtureFraction& flamelet
)
{
    Info<< "     c-o Solve laplace equation to get initial solution\n";

    //- Solve the Laplace Equation to get initial (linear) profile
    //  We only need the boundary species
    
    //- Oxidizer species 
    const wordList& speciesO = flamelet.speciesOxidizer();

    //- Fuel species
    const wordList& speciesF = flamelet.speciesFuel();

    //- Species word list
    wordList species_ = speciesO;

    //- Extend with fuel species
    forAll(speciesF, s)
    {
        species_.push_back(s);
    }

    //- Calculate laplace equation with diffusion coefficient D = 1
    //  Note: Fourier No. < 0.5 -> Fo = dt * D / dZ^2
    //  We have uniform spaced dZ till now
    const scalar dZ = scalar(1) / (flamelet.nZPoints() -1);

    const scalar dt = scalar(0.49) * pow(dZ, 2);

    const int nZ = flamelet.nZPoints();

    const scalarField& Zvalues = flamelet.Z();

    //- Old field
    scalarField Yold(nZ, 0);

    //- New Y field
    scalarField Ynew(nZ, 0);

    forAll(species_, s)
    {
        size_t nIter{0};

        //- FinalResidual
        scalar finalResidual{0};

        do
        {
            //- Iteration
            ++nIter;

            //- Store old iteration Y field
            Yold = flamelet.Y(s);

            //- Copy boundary conditions
            Ynew[0] = Yold[0];
            Ynew[nZ-1] = Yold[nZ-1];

            //- Solve laplace equation for each discrete point
            for (int i = 1; i < nZ-1; ++i)
            {
                Ynew[i] =
                    (
                        FDMLapacian2ndOrder
                        (
                            Yold[i-1],
                            Yold[i],
                            Yold[i+1],
                            Zvalues[i-1],
                            Zvalues[i],
                            Zvalues[i+1]
                        )
                    ) * dt + Yold[i];
            }

            //- Update the mass fraction
            flamelet.updateY(s, Ynew);

            finalResidual = residual(Yold, Ynew);

            //- Exit based on mass conservation TODO
        }
        while (finalResidual > 1e-9);

        //- Update the other fields
        flamelet.updateFields();

        Info<< "         ++ Solved for species " << s << ", Converged after "
            << nIter << ", Residual = " << finalResidual << "\n";
    }

}


template<class Type>
void AFC::Numerics<Type>::solveAdiabaticFlamelet
(
    MixtureFraction& flamelet
)
{
    Info<< "     c-o Solve the flamelet equation for the adiabatic flamelet\n";

    // --------- NEW 
    //
    //
    //
    //


    //- Iteration counter
    size_t nIter{0};

    scalar tEnd = 10;
    scalar deltaTDiff = 0.05;
    scalar deltaTChem = 1e-8;
    scalar dt{0.5};

    //- Test case
    //  Use one point and use it as ideal reactor
    const map<word, scalar>& c = flamelet.C(2);
    const wordList& species = flamelet.species();
    const scalar T = flamelet.T(2);

    ode_->solve(T, species, c);

    //- Solve the flamelet equation. Therefore we split the calculation into
    //  chemistry solving with the chemistry delta till we reach the diffusion
    //  delta. After that we go on in time for the next time step
    //
    //  Time:
    //      tEnd := Simulation end time 
    //      deltaTDiff := Diffusion time step (for flamelet equation)
    //      deltaTChem := Chemistry time step (sub time step for chemistry)
    //    
    //  Chemistry is solved using the Seulex algorithm
    //
    //  
    //    

    //- TH::Solver chemistry.solve()
    //  {
    //- Run time (solve till we reach end time)
    
    /*while (deltaTDiff < tEnd)
    {
        //- Store old chemistry time step
        scalar deltaTChem0 = deltaTChem;

        //- Solve the chemistry and return the maximum time step for the 
        //  chemistry during this calculation (used for the next iteration)
        deltaTChem = solveChemistry(dt, flamelet);


        //-
    }
    //- TH::Solver chemistry.solve()
    //  }
    //  */
}


template<class Type>
AFC::scalar AFC::Numerics<Type>::solveChemistry
(
    const scalar dt,
    MixtureFraction& flamelet
)
{

    //- TH::ChemistryModel
    //{
    scalar deltaTMin{1e6};

    //- Temperature and pressure
    const scalarField& T = flamelet.T();
    const scalar p = flamelet.p();

    //- Copy of concentration field
    const List<map<word, scalar> > c0 = flamelet.C();
    
    //- The density field [g/m^3]
    const scalarField& rho = flamelet.rho();

    //- Solve the chemistry for each discrete mixture fraction point
    forEach(T, Zi)
    {
        //- Skip the first and last point (pure oxidizer and fuel)
        if ((Zi > 0) && (Zi < T.size()))
        {
            //- Get values at the discrete point (pressure is constant)
            scalar Ti = T[Zi];
            scalar rhoi = rho[Zi];
            const map<word, scalar>& c0i = c0[Zi];

            //- New concentration field that is been updated
            map<word, scalar> ci = c0[Zi];

            //ode_->solve(Ti, p, ci, dt, deltaTChem_[Zi]);

        }
    }

    return deltaTMin;
    //- TH::ChemistryModel
    //}
}


template<class Type>
AFC::scalar AFC::Numerics<Type>::FDMLapacian2ndOrder
(
    const scalar& phi_l,
    const scalar& phi,
    const scalar& phi_r,
    const scalar& Z_l,
    const scalar& Z,
    const scalar& Z_r
) const
{
    //- Return the solution of the laplacian term
    //  l == left position
    //  r == right position
    //  no subscript is the value in the center
    //  Derivation given in [Holzmann]
    return
    (
        //- Enumerator
        (
            (phi_r - phi) / (Z_r - Z)
          - (phi - phi_l) / (Z - Z_l)
        )
        //- Denominator
      / (scalar(0.5) * (Z_r - Z_l))
    );
}


template<class Type>
void AFC::Numerics<Type>::solveFlamelet
(
    MixtureFraction& flamelet,
    const scalar& chi,
    const scalar& dt
)
{
    //- Save old fields
    const List<map<word, scalar> >& Yold = flamelet.Y();
    const List<map<word, scalar> >& Xold = flamelet.X();
    const List<map<word, scalar> >& Cold = flamelet.C();
    const scalarField& Told = flamelet.T();

    //- Space discretization
    const scalarField& Zvalues = flamelet.Z();
    const int nZ = flamelet.nZPoints();

    //- Species that need to be solved
    const wordList species_ = flamelet.species();

    //- Solve species equation for mass fraction
    //  Time derivation => 1nd Order Euler
    //  Laplace derivation => 2nd Order Linear for non-equal spaced distances
    //  Source term omega handled explicit (for the first implementation)
    //  Equation is solved explicit (implicit should be implemented)

    //- Diffusion coefficient (chi / 2)
    const scalar D = chi/scalar(2);

    //- Source term omega
    const scalar omega = 0;

    //- Unset updated
    flamelet.updatedFields(false);
    
    //- Solve the flamelet equation for each species
    //  Derivation is given in [Holzmann]

    //- a) Calculate and update the fields
    for (int i = 1; i < nZ-1; ++i)
    {
        //- 1) Mean molecular weight
        flamelet.updateMeanMW();

        //- 2) Density
        flamelet.updateRho();
    }

    forAll(species_, s)
    {
        //- Temporar field that stores the new values
        scalarField Ynew(nZ, 0);

        //- Copy boundary conditions
        Ynew[0] = Yold[0].at(s);
        Ynew[nZ-1] = Yold[nZ-1].at(s);

        //- Solve for each discrete point (first and last points are BC)
        for (int i = 1; i < nZ-1; ++i)
        {
            Ynew[i] =
                (
                    D
                  * FDMLapacian2ndOrder
                    (
                        Yold[i-1].at(s),
                        Yold[i].at(s),
                        Yold[i+1].at(s),
                        Zvalues[i-1],
                        Zvalues[i],
                        Zvalues[i+1]
                    )
                  + omega
                ) * dt + Yold[i].at(s);
        }

        //- Update the fields (could be done with references faster)
        flamelet.updateY(s, Ynew);
        flamelet.updateX();
        flamelet.updateC(); 
    }

    //- Solve the flamelet equation for the temperature
    {
        //- Temporar field that stores the new temperatures 
        scalarField Tnew(nZ, 0);

        //- Copy boundary conditions
        Tnew[0] = Told[0];
        Tnew[nZ-1] = Told[nZ-1];

        //- Further source terms

        //- Solve the temperature equation
        for (int i = 1; i < nZ-1; ++i)
        {
            Tnew[i] =
                (
                    D
                  * FDMLapacian2ndOrder
                    (
                        Told[i-1],
                        Told[i],
                        Told[i+1],
                        Zvalues[i-1],
                        Zvalues[i],
                        Zvalues[i+1]
                    )
                  + omega
                ) * dt + Told[i];
        }

        //- Update the field
        flamelet.updateT(Tnew);
    }
}


/*
void AFC::Numerics<Type>::jacobian
(
    MixtureFraction& mf
)
{
    if (debug)
    {
        Info<< " --> AFC::Numerics<Type>::jacobian" << endl;
    }
}
*/

/*void calculate
(
    lookUpTable& lut,
    const scalar& sDR,
    const scalar& defect,
    const unsigned int nDisPoints
)
{
    //- TODO
    //- Calculate reaction rate factors k for each reaction

    for (unsigned int point=0; point <= nDisPoints; point++)
    {
//        Info<< "Z = " << point << "\n";
        //- Object of discrete mixture fraction
        MixtureFraction& dMF = lut[defect][sDR][point];

        //- Temperature at discrete point Z
        const scalar& T = dMF.T();

        //- Mol fractions at discrete point
//        const map<word, scalar>& speciesMol = dMF.mol();
        {
            //- a) calculate mean molecular weight MW (using mol)
            //dMF.calculateMeanMW("mol");

            //- b) calculate mean heat capacity cp [J/mol/K]
            //dMF.calculateMeanCp(T);

            //- c) calculate mean enthalpy H [J/mol]
            //dMF.calculateMeanH(T);

            //- d) calculate mean entropy S [J/mol/K]
            //dMF.calculateMeanS(T);

            //- e) calculate mean free gibbs energy
            {
             //   const scalar& H = dMF.H();

             //   const scalar& S = dMF.S();

             //   dMF.calculateMeanG(H, S, T);
            }
        }
    }
}

        //- Calculate the source term of species omega
        {
            //- Actual concentration of species at point Zi
            //const map<word, scalar>& con1 = dMF.con();

            //- Copy of species concentration
            //map<word, scalar> con = con1;

            //dMF.calculateOmega(T, con);

            //- Update kf and kb using the new T field

            //- Calc for each species the source term omega
            //chem_.omega(speciesMol);
            //
            //I
        }
    }
}*/


template<class Type>
AFC::scalar AFC::Numerics<Type>::residual
(
    const scalarField& oldField,
    const scalarField& newField
) const
{
    scalarField residualField(newField.size(), 0);

    for (unsigned int i = 0; i < newField.size(); ++i)
    {
        residualField[i] = fabs(newField[i] - oldField[i]);
    }

    return max(residualField);
}


template<class Type>
AFC::scalar AFC::Numerics<Type>::max
(
    const scalarField& sF
) const
{
    scalar maxValue{0};

    forAll(sF, value)
    {
        if (value >= maxValue)
        {
            maxValue = value;
        }
    }

    return maxValue;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template class AFC::Numerics<AFC::Rosenbrock>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
/*
                //- e) if species is in reaction we go on in calculating the
                //  derivative
                if (cont)
                {
                    //- f) Product and educt species
                    const wordList& prodSpecies = flamelet.speciesProducts(r);
                    const wordList& educSpecies = flamelet.speciesEducts(r);

                    //- g) Concentrations of species at discrete point
                    const map<word, scalar>& con = flamelet.C(Z);

                    //- Where is the derivative species?
                    bool inEducSite{false};
                    bool inProdSite{false};

                    forAll(prodSpecies, s)
                    {
                        if (s == species2)
                        {
                            inProdSite = true;
                            break;
                        }
                    }

                    forAll(educSpecies, s)
                    {
                        if (s == species2)
                        {
                            inEducSite = true;
                            break;
                        }
                    }

                    Info<< "In product site: " << inProdSite << "\n";
                    Info<< "In educt   site: " << inEducSite << "\n";

                    //- h) Derivative species is in reaction as educt
                    if (inEducSite)
                    {
                        scalar educ{0};

                        //- i) Get the stochiometric coeffs
                        const map<word, int>& nuProd = flamelet.nuProducts(r);
                        const map<word, int>& nuEduc = flamelet.nuEducts(r);

                        //- j) Stochiometric factor of the derivative species
                        const scalar nuDS = nuEduc.at(species2);

                        Info<< "  The derivative species is of order " << nuDS << "\n";
                        Info<< "  Educt = ";

                        //- k) nuDS == 1 (linear) nuDS == 2 (non-linear)
                        if (fabs(nuDS) == 1)
                        {
                            forAll(educSpecies, s)
                            {
                                if (s != species2)
                                {
                                    Info<< "["<< s << "]^" << fabs(nuEduc.at(s)) << " * "; 
                                    educ *= pow(con.at(s), fabs(nuEduc.at(s)));
                                }
                            }

                        }
                        else if (fabs(nuDS) == 2)
                        {
                            //- Avoid repeating
                            bool repeating{false};

                            forAll(educSpecies, s)
                            {
                                if (s != species2)
                                {
                                    Info<< "["<< s << "]^" << fabs(nuEduc.at(s)) << " * "; 
                                    educ *= pow(con.at(s), fabs(nuEduc.at(s)));
                                }
                                else
                                {
                                    if (!repeating)
                                    {
                                    Info<< fabs(nuEduc.at(s)) << "["<< s << "]^"
                                        << fabs(nuEduc.at(s))-1 << " * "; 

                                    //- First derivative of the term
                                    educ *= (fabs(nuEduc.at(s)))
                                        * pow(con.at(s), fabs(nuEduc.at(s))-1);
                                    repeating = true;
                                    }
                                }
                            }
                        }
                        else
                        {
                            ErrorMsg
                            (
                                "The stochiometric factor is not defined\n",
                                __FILE__,
                                __LINE__
                            );
                        }
                        Info<<"\n";

                        //- l) get forward reaction rate kf
                        const scalar kf = flamelet.kf(r, Z);

                        //- Get pre-factor nu'' - nu' :: based on the fact that we already
                        //  know the right value, we just have to check if the species
                        //  is within the product or educt side
                        scalar nuSpecies{0};

                        if (nuEduc.count(species2))
                        {
                            nuSpecies = nuEduc.at(species2);
                        }
                        else
                        {
                            nuSpecies = nuProd.at(species2);
                        }

                        //- n) Add the value to the Jacobian element
                        Info<< "---> " << nuSpecies << " * " << kf << " * " << educ 
                            << " = " << nuSpecies * -1 * kf * educ <<  "\n";
                        Jij += nuSpecies * kf * educ;
                    }

                    //- h) Derivative species is in reaction as product
                    if (inProdSite)
                    {
                        scalar prod{0};

                        //- i) Get the stochiometric coeffs
                        const map<word, int>& nuProd = flamelet.nuProducts(r);
                        const map<word, int>& nuEduc = flamelet.nuEducts(r);

                        //- j) Stochiometric factor of the derivative species
                        const scalar nuDS = nuProd.at(species2);

                        Info<< "  The derivative species is of order " << nuDS << "\n";
                        Info<< "  Prod = ";

                        //- k) nuDS == 1 (linear) nuDS == 2 (non-linear)
                        if (fabs(nuDS) == 1)
                        {
                            forAll(prodSpecies, s)
                            {
                                if (s != species2)
                                {
                                    Info<< "["<< s << "]^" << nuProd.at(s) << " * "; 
                                    prod *= pow(con.at(s), nuProd.at(s));
                                }
                            }

                        }
                        else if (fabs(nuDS) == 2)
                        {
                            //- Avoid repeating
                            bool repeating{false};

                            forAll(prodSpecies, s)
                            {
                                if (s != species2)
                                {
                                    Info<< "["<< s << "]^" << nuProd.at(s) << " * "; 
                                    prod *= pow(con.at(s), nuProd.at(s));
                                }
                                else
                                {
                                    if (!repeating)
                                    {
                                    //- First derivative of the term
                                    Info<< nuProd.at(s) << "["<< s << "]^" << nuProd.at(s)-1 << " * "; 
                                    prod *= (nuProd.at(s) - 1)
                                        * pow(con.at(s), nuProd.at(s)-1);
                                    repeating = true;
                                    } 
                                }
                            }
                        }
                        else
                        {
                            ErrorMsg
                            (
                                "The stochiometric factor is not defined\n",
                                __FILE__,
                                __LINE__
                            );
                        }
                        Info<<"\n";

                        //- l) get backward reaction rate kb
                        const scalar kb = flamelet.kb(r, Z);

                        //- Get pre-factor nu'' - nu' :: based on the fact that we already
                        //  know the right value, we just have to check if the species
                        //  is within the product or educt side
                        scalar nuSpecies{0};

                        if (inEducSite)
                        {
                            nuSpecies = nuEduc.at(species2);
                        }
                        else if (inProdSite) 
                        {
                            nuSpecies = nuProd.at(species2);
                        }
                        else if (inEducSite && inProdSite)
                        {
                            nuSpecies = nuProd.at(species2) + nuProd.at(species2);
                        }

                        //- n) Add the value to the Jacobian element
                        Info<< "---> " << nuSpecies << " * - " << kb << " * " << prod
                            << " = " << nuSpecies * -1 * kb * prod <<  "\n";
                        Jij += nuSpecies * -1 * kb * prod; 
                    }
*/
