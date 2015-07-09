/*---------------------------------------------------------------------------*\
| =========                |                                                  |
| \\      /  A utomatic    | Holzmann-cfd                                     |
|  \\    /   F lamelet     | Version: 1.0                                     |
|   \\  /    C onstructor  | Web: www.Holzmann-cfd.de                         |
|    \\/                   |                                                  |
\*---------------------------------------------------------------------------*/
/*
»
»
»
»
»
»
»
»
»
»
\*---------------------------------------------------------------------------*/
//- system headers


//- user def. headers
#include "functions.hpp"
#include "../transport/transport.hpp"

//- non class functions

const unsigned int mixtureFractionPoints
(
    const normalString& fileName
)
{
    //- output
    std::cout << "Reading mixtureFractionPoints (" << fileName << ")\n";

    //- read the whole file and store it
    stringField fileContent = Transport::openFile(fileName);

    //- keyword found
    bool found{false};

    //- loop through the fileContent
    forAll(fileContent, line)
    {
        //- STEP 1: remove all whitespaces
        stringField tmp = Transport::splitString(fileContent[line]);

        //- if not empty, and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            if (tmp[0] == "mixtureFractionPoints")
            {
                std::cout << " ++ Mixture fraction points: " << tmp[1] << "\n";
                return stoi(tmp[1]);
            }
        }
    }

    //- if keyword not found abort
    if (!found)
    {
        std::cerr<< "\n ++ ERROR: Keyword \"mixtureFractionPoints\" not "
                 << "found in file " << fileName << "."
                 << __FILE__ << " line " << __LINE__ << std::endl;
        std::terminate();
    }

    //- dummy
    return -1;
}


//- scalar dissipation rates
const scalarField scalarDissipationRates
(
    const normalString& fileName
)
{
    //- output
    std::cout << "Reading scalar dissipation rates [Hz] (" << fileName << ")\n";

    //- read the whole file and store it
    stringField fileContent = Transport::openFile(fileName);

    bool scalarDissipationRateFound{false};
    bool inDictionary{false};

    scalarField scalarDissipationRate;

    //- loop through the fileContent
    forAll(fileContent, line)
    {
        //- STEP 1: remove all whitespaces
        stringField tmp = Transport::splitString(fileContent[line]);

        //- if not empty, and no comment, proceed
        if
        (
            !tmp.empty()
         && tmp[0][0] != '!'
        )
        {
            //- moleFractionFuel and moleFractionOxidizer
            if
            (
                tmp[0] == "scalarDissipationRates"
            )
            {
                scalarDissipationRateFound = true;
            }
            //- scalar dissipation rate dictionary
            else if
            (
                scalarDissipationRateFound
             && !inDictionary
            )
            {
                if (tmp[0] == "{")
                {
                    inDictionary = true;
                }
            }
            else if
            (
                scalarDissipationRateFound
             && inDictionary
            )
            {
                if (tmp[0] == "}")
                {
                    forAll(scalarDissipationRate, i)
                    {
                        std::cout << "  ++ " << scalarDissipationRate[i] << " Hz\n";
                    }
                    std::cout << "\n";

                    return scalarDissipationRate;
                }
                else
                {
                    //- if no number this line produces an error
                    scalarDissipationRate.push_back(stod(tmp[0]));
                }
            }
        }
    }

    //- dummy
    return std::vector<double>(0);
}


//        //- set mass fraction for species
//        //  could be done in a better way
//        scalar MMW_fuel{0};
//        scalar MMW_oxidizer{0};
//
//        //- calculate mean molecular weight of fuel and oxidizer
//        forAll(species, id)
//        {
//            if (species[id].fuel())
//            {
//                MMW_fuel += species[id].MW() * species[id].Xf();
//            }
//            if (species[id].oxidizer())
//            {
//                MMW_oxidizer += species[id].MW() * species[id].Xo();
//            }
//        }
//
//        //- calculate mass fraction Y
//        forAll(species, id)
//        {
//            if (species[id].fuel())
//            {
//                species[id].setYf
//                (
//                    species[id].MW() /
//                    MMW_fuel *
//                    species[id].Xf()
//                );
//            }
//            if (species[id].oxidizer())
//            {
//                species[id].setYo
//                (
//                    species[id].MW() /
//                    MMW_oxidizer *
//                    species[id].Xo()
//                );
//            }
//        }
//
//        //- check if X = 1  and Y = 1
//        scalar X{0};
//        scalar Y{0};
//        std::cout<< "Checking chemical composition of fuel:\n"
//                 << "------------------------------------------------------\n";
//        forAll(species, id)
//        {
//            if (species[id].fuel())
//            {
//                std::cout<< " ++ " << species[id].name()
//                         << "\t\t X = " << species[id].Xf()
//                         << "\t Y = " << species[id].Yf()
//                         << "\t MW = " << species[id].MW() << "\n";
//                X += species[id].Xf();
//                Y += species[id].Yf();
//            }
//        }
//        std::cout<< " ++ Fuel mol fraction = " << X << "\n";
//        std::cout<< " ++ Fuel mass fraction = " << Y << "\n\n";
//
//        //- if sum <> 1
//        if (X != 1. || Y != 1.)
//        {
//            std::cerr<< "\n ++ " << "Fuel mol/mass fraction "
//                     << "is not equal to 1"
//                     << "\n ++ Error occur in file " << __FILE__
//                     << " line " << __LINE__ << std::endl;
//            std::terminate();
//        }
//
//        //- reset X and Y
//        X = 0;
//        Y = 0;
//        std::cout<< "\nChecking chemical composition of oxidizer:\n"
//                 << "------------------------------------------------------\n";
//        forAll(species, id)
//        {
//            if (species[id].oxidizer())
//            {
//                std::cout<< " ++ " << species[id].name()
//                         << "\t\t X = " << species[id].Xo()
//                         << "\t Y = " << species[id].Yo()
//                         << "\t MW = " << species[id].MW() << "\n";
//                X += species[id].Xo();
//                Y += species[id].Yo();
//            }
//        }
//        std::cout<< " ++ Oxidizer mol fraction = " << X << "\n";
//        std::cout<< " ++ Oxidizer mass fraction = " << Y << "\n\n\n";
//
//        //- if sum <> 1
//        if (X != 1. || Y != 1.)
//        {
//            std::cerr<< "\n ++ " << "Oxidizer mol/mass fraction "
//                     << "is not equal to 1"
//                     << "\n ++ Error occur in file " << __FILE__
//                     << " line " << __LINE__ << std::endl;
//            std::terminate();
//        }
//    }
//}
//

//    //- open file and return the whole file as string array
//    const stringField openFile(const normalString& fileName)
//    {
//        std::ifstream file;
//
//        //- open file
//        file.open(fileName.c_str(), std::ios::in);
//
//        //- check file
//        if(!file.good())
//        {
//            std::cerr<< "\n ++ ERROR: Could not open file \""
//                     << fileName << "\"..."
//                     << "\n ++ Error occur in file " << __FILE__
//                     << " line " << __LINE__ << std::endl;
//            std::terminate();
//        }
//
//        //- stored line
//        normalString fileLine;
//
//        //- file content
//        stringField fileContent(0);
//
//        //- read the whole file
//        while(!file.eof())
//        {
//            std::getline(file, fileLine);
//            fileContent.push_back(fileLine);
//        }
//        return fileContent;
//    }
//
//
//    //- calculate stoichiometric mixture fraction Zst
//    scalar stochiometricMF
//    (
//        const std::vector<Species>& species
//    )
//    {
//        std::cout<< "Calculate stochiometric mixture fraction Zst:\n"
//                 << "------------------------------------------------------\n";
//
//        //- mean mol fraction in [g/mol]
//        scalar MMW_fuel{0};
//        scalar MMW_oxidizer{0};
//
//        //- Step 1
//        //  Calculate mean molecular weight of fuel and oxidizer
//        forAll(species, id)
//        {
//            if (species[id].fuel())
//            {
//                MMW_fuel += species[id].MW() * species[id].Xf();
//            }
//            if (species[id].oxidizer())
//            {
//                MMW_oxidizer += species[id].MW() * species[id].Xo();
//            }
//        }
//
//        std::cout<< "  ++ Mean molecular weight fuel: " << MMW_fuel << "\n"
//                 << "  ++ Mean molecular weight oxidizer: " << MMW_oxidizer << "\n\n";
//
//        //- element mass fraction in fuel stream
//        scalar Zc_F{0};
//        scalar Zh_F{0};
//        scalar Zo_F{0};
//
//        //- element mass fraction in oxidizer stream
//        scalar Zc_O{0};
//        scalar Zo_O{0};
//        scalar Zh_O{0};
//        scalar Zn_O{0};
//
//
//        //- Step 2
//        //  Calculate element mass fraction Zc, Zh and Zo in fuel
//        //  Z_i = M_i / MMW_xxx * sum ( a_ij * X_j )
//        //      ++ Z_i:  element mass fraction of element i [-]
//        //      ++ M_i:  molecular weight of element i [g/mol]
//        //      ++ a_ij: amount of element i in species j [-]
//        //      ++ X_j:  mol fraction of species j [-]
//
//        std::cout<< "\nComposition \tC\tH\tO\tN\n"
//                 << "------------------------------------------------------\n";
//        forAll(species, id)
//        {
//            //- TODO use inert and placeholder
//            scalar nCAtoms_F{0};
//            scalar nHAtoms_F{0};
//            scalar nOAtoms_F{0};
//            scalar nNAtoms_F{0};    // INERT not used
//
//            scalar nOAtoms_O{0};
//            scalar nNAtoms_O{0};    // INERT
//            scalar nCAtoms_O{0};    // placeholder
//            scalar nHAtoms_O{0};    // H in oxidizer
//
//            normalString aP, nP;
//
//            if (species[id].fuel() || species[id].oxidizer())
//            {
//                //- TODO (two character elements)
//                normalString tmp = species[id].name();
//
//                for (unsigned int i=0; i<tmp.length(); i++)
//                {
//                    //- actual position
//                    aP = tmp[i];
//
//                    //- next position
//                    if (i+1 <= tmp.length())
//                    {
//                        nP = tmp[i+1];
//                    }
//                    else
//                    {
//                        nP = aP;
//                    }
//                    //- only fuel
//                    if (species[id].fuel())
//                    {
//                        //- C atom
//                        if (aP == "C")
//                        {
//                            //- C is last character, no number
//                            //  one C atom
//                            if (aP == nP)
//                            {
//                                nCAtoms_F += 1;
//                            }
//                            //- after C there is a character
//                            else
//                            {
//                                //- if next char is a number its the amount of
//                                //  atoms of C
//                                if
//                                (
//                                    nP.find_first_of("0123456789") !=
//                                    std::string::npos
//                                )
//                                {
//                                    nCAtoms_F += stod(nP);
//                                }
//                                //- if next char is a character the atom is
//                                //  only used once
//                                else
//                                {
//                                    nCAtoms_F++;
//                                }
//                            }
//                        }
//                        //- H atom
//                        if (aP == "H")
//                        {
//                            //- H is last character, no number
//                            //  one H atom
//                            if (aP == nP)
//                            {
//                                nHAtoms_F++;
//                            }
//                            //- after H there is a character
//                            else
//                            {
//                                //- if next char is a number its the amount of
//                                //  atoms of H
//                                if
//                                (
//                                    nP.find_first_of("0123456789") !=
//                                    std::string::npos
//                                )
//                                {
//                                    nHAtoms_F += stod(nP);
//                                }
//                                //- if next char is a character the atom is
//                                //  only used once
//                                else
//                                {
//                                    nHAtoms_F++;
//                                }
//                            }
//                        }
//                        //- O atom
//                        if (aP == "O")
//                        {
//                            //- O is last character, no number
//                            //  one O atom
//                            if (aP == nP)
//                            {
//                                nOAtoms_F++;
//                            }
//                            //- after O there is a character
//                            else
//                            {
//                                //- if next char is a number its the amount of
//                                //  atoms of O
//                                if
//                                (
//                                    nP.find_first_of("0123456789") !=
//                                    std::string::npos
//                                )
//                                {
//                                    nOAtoms_F += stod(nP);
//                                }
//                                //- if next char is a character the atom is
//                                //  only used once
//                                else
//                                {
//                                    nOAtoms_F++;
//                                }
//                            }
//                        }
//                    }
//                    if (species[id].oxidizer())
//                    {
//                        //- C atom
//                        if (aP == "C")
//                        {
//                            //- C is last character, no number
//                            //  one C atom
//                            if (aP == nP)
//                            {
//                                nCAtoms_O += 1;
//                            }
//                            //- after C there is a character
//                            else
//                            {
//                                //- if next char is a number its the amount of
//                                //  atoms of C
//                                if
//                                (
//                                    nP.find_first_of("0123456789") !=
//                                    std::string::npos
//                                )
//                                {
//                                    nCAtoms_O += stod(nP);
//                                }
//                                //- if next char is a character the atom is
//                                //  only used once
//                                else
//                                {
//                                    nCAtoms_O++;
//                                }
//                            }
//                        }
//                        //- O atom
//                        if (aP == "O")
//                        {
//                            //- O is last character, no number
//                            //  one O atom
//                            if (aP == nP)
//                            {
//                                nOAtoms_O++;
//                            }
//                            //- after O there is a character
//                            else
//                            {
//                                //- if next char is a number its the amount of
//                                //  atoms of O
//                                if
//                                (
//                                    nP.find_first_of("0123456789") !=
//                                    std::string::npos
//                                )
//                                {
//                                    nOAtoms_O += stod(nP);
//                                }
//                                //- if next char is a character the atom is
//                                //  only used once
//                                else
//                                {
//                                    nOAtoms_O++;
//                                }
//                            }
//                        }
//                        //- H atom
//                        if (aP == "H")
//                        {
//                            //- H is last character, no number
//                            //  one H atom
//                            if (aP == nP)
//                            {
//                                nHAtoms_O++;
//                            }
//                            //- after H there is a character
//                            else
//                            {
//                                //- if next char is a number its the amount of
//                                //  atoms of H
//                                if
//                                (
//                                    nP.find_first_of("0123456789") !=
//                                    std::string::npos
//                                )
//                                {
//                                    nHAtoms_O += stod(nP);
//                                }
//                                //- if next char is a character the atom is
//                                //  only used once
//                                else
//                                {
//                                    nHAtoms_O++;
//                                }
//                            }
//                        }
//                        //- N atom
//                        if (aP == "N")
//                        {
//                            //- N is last character, no number
//                            //  one N atom
//                            if (aP == nP)
//                            {
//                                nNAtoms_O++;
//                            }
//                            //- after N there is a character
//                            else
//                            {
//                                //- if next char is a number its the amount of
//                                //  atoms of N
//                                if
//                                (
//                                    nP.find_first_of("0123456789") !=
//                                    std::string::npos
//                                )
//                                {
//                                    nNAtoms_O += stod(nP);
//                                }
//                                //- if next char is a character the atom is
//                                //  only used once
//                                else
//                                {
//                                    nNAtoms_O++;
//                                }
//                            }
//                        }
//                    }
//
//                }
//                //- output fuel composition
//                if (species[id].fuel())
//                {
//                    std::cout<< tmp << "\t\t" << nCAtoms_F << "\t"
//                             << nHAtoms_F << "\t" << nOAtoms_F << "\t"
//                             << nNAtoms_F << "\t (F)\n";
//                }
//
//                //- output oxidzier composition
//                if (species[id].oxidizer())
//                {
//                    std::cout<< tmp << "\t\t" << nCAtoms_O << "\t"
//                             << nHAtoms_O << "\t" << nOAtoms_O << "\t"
//                             << nNAtoms_O << "\t (O)\n";
//                }
//
//                //- element mass fraction (without M_i/M)
//                Zc_F += nCAtoms_F * species[id].Xf();
//                Zh_F += nHAtoms_F * species[id].Xf();
//                Zo_F += nOAtoms_F * species[id].Xf();
//                Zo_O += nOAtoms_O * species[id].Xo();
//                Zc_O += nCAtoms_O * species[id].Xo();
//                Zh_O += nHAtoms_O * species[id].Xo();
//                Zn_O += nNAtoms_O * species[id].Xo();
//            }
//        }
//
//        //- access to database
//        Elements element;
//
//        Zc_F *= element.atomicWeight("C") / MMW_fuel;
//        Zh_F *= element.atomicWeight("H") / MMW_fuel;
//        Zo_F *= element.atomicWeight("O") / MMW_fuel;
//
//        Zc_O *= element.atomicWeight("C") / MMW_oxidizer;
//        Zh_O *= element.atomicWeight("H") / MMW_oxidizer;
//        Zo_O *= element.atomicWeight("O") / MMW_oxidizer;
//        Zn_O *= element.atomicWeight("N") / MMW_oxidizer;
//
//        std::cout<< "\n\nElement mass fraction:\n"
//                 << "------------------------------------------------------\n"
//                 << "  ++ Zc\t\t " << Zc_F << "\t\t (F)\n"
//                 << "  ++ Zh\t\t " << Zh_F << "\t\t (F)\n"
//                 << "  ++ Zo\t\t " << Zo_F << "\t\t (F)\n"
//                 << "  ++ Sum of element mass fraction: " << Zc_F + Zh_F + Zo_F
//                 << "\t (F)\n\n"
//                 << "  ++ Zc\t\t " << Zc_O << "\t\t (O)\n"
//                 << "  ++ Zh\t\t " << Zh_O << "\t\t (O)\n"
//                 << "  ++ Zo\t\t " << Zo_O << "\t\t (O)\n"
//                 << "  ++ Zn\t\t " << Zn_O << "\t\t (O)\n"
//                 << "  ++ Sum of element mass fraction: "
//                 << Zc_O + Zh_O + Zo_O + Zn_O
//                 << "\n\n";
//
//        //- step 3
//        //  calculate omin
//        scalar omin{0};
//        scalar factor_Zc{0};
//        scalar factor_Zh{0};
//
//        //- N. Peters
//        factor_Zc = element.atomicWeight("O")*2/element.atomicWeight("C");
//        factor_Zh = element.atomicWeight("O")*2/(element.atomicWeight("H")*4);
//
//        //- omin
//        //  [kgO2/kgF]
//        //  [gO2/gF]
//        //  TODO
//        //  >>> if H or C in oxidizer we need more O2 <<< ???
//        //  >>> also not included in libOpenSMOKE??? <<<
//        omin = factor_Zc * Zc_F + factor_Zh * Zh_F - Zo_F;
//        std::cout<< "  ++ omin: " << omin << "\n";
//
//        //- step 4
//        //- calculate the amount of moles O per
//        //  mole F at stochiometric condition
//        //  o_st [molO/molF]
//        scalar o_st = omin*MMW_fuel / (element.atomicWeight("O")*2) / Zo_O;
//
//        //- step 5
//        //  calculate mean moleculare weight of fuel and oxidizer
//        //  at stochiometric conditions
//        scalar MMW_st = MMW_fuel + o_st * element.atomicWeight("O")*2;
//
//        //- step 6
//        //  calculate elementar mass fraction Zo at stochiometric condition
//        scalar Zo_st = o_st * element.atomicWeight("O")*2 / MMW_st;
//
//        std::cout<< "  ++ Stochiometric mixture fraction Zst: " << Zo_st << "\n";
//
//        return Zo_st;
//    }
//
//    //- calculate adiabatic enthalpy of fuel and oxidizer
//    //  + 0 mean fuel
//    //  + 1 mean oxidizer
//    //  [J/kg]
//    scalar adiabaticEnthalpy
//    (
//        const std::vector<Species>& species,
//        const int i
//    )
//    {
//        scalar hf_a{0};
//        scalar ho_a{0};
//
//        forAll(species, id)
//        {
//            // hf_a     [J/kg]
//            // h()      [J/mol]
//            // MW()     [g/mol]
//            // Xf()     [-]
//            // 1000     g/kg
//            if (species[id].fuel())
//            {
//                hf_a += species[id].h(species[id].Tf()) *
//                        species[id].Xf() /
//                        species[id].MW() *
//                        1000;
//            }
//            if (species[id].oxidizer())
//            {
//                ho_a += species[id].h(species[id].To()) *
//                        species[id].Xo() /
//                        species[id].MW() *
//                        1000;
//            }
//        }
//
//        //- fuel
//        if (i == 0)
//        {
//            return hf_a;
//        }
//        //- oxidizer
//        else if (i == 1)
//        {
//            return ho_a;
//        }
//        else
//        {
//            std::cerr<< "\n ++ ERROR: wrong value for quantifier <F=0/O=1>"
//                     << "\n ++ Error occur in file " << __FILE__
//                     << " line " << __LINE__ << std::endl;
//            std::terminate();
//        }
//    }
//
//    //- calculate adiabatic flame temperature [K]
//    //  for stochiometric conditions
//    /*scalar adiabateFlameTemperature
//    (
//        const scalar& Zst,
//        const std::vector<Species>& species
//    )
//    {
//        scalar hu0{0};
//        int idH2O{-1};
//        int idN2{-1};
//
//        forAll(species, id)
//        {
//            if (species[id].name() == "O2")
//            {
//                idH2O = id;
//            }
//            if (species[id].name() == "N2")
//            {
//                idN2 = id;
//            }
//        }
//
//        std::cout << "---> " << species[idH2O].h(298)/1000 << "\n";
//
//        //- ref temperature
//        scalar Tref{298};
//
//        //- point 1
//        scalar T1{500};
//        scalar Tm1{0};
//        scalar hb1{0};
//        scalar cpm1{0};
//
//        //- point 2
//        scalar T2{5000};
//        scalar Tm2{0};
//        scalar hb2{0};
//        scalar cpm2{0};
//
//        //- mid point
//        scalar T3{0};
//        scalar Tm3{0};
//        scalar hb3{0};
//        scalar cpm3{0};
//
//        //- normalized values
//        scalar normalized1{0};
//        scalar normalized2{0};
//        scalar normalized3{0};
//
//        scalar href=545170.87;
//
//        int i{0};
//        scalar T{298};
//        do
//        { break;
//            //- arithmetic mean of T1 and T2
//            //  based on T=298K
//            Tm1 = (T1+Tref)/2;
//            Tm2 = (T2+Tref)/2;
//
//            //- mid point temperature
//            T3 = (T1+T2)/2;
//            Tm3 = (T3+Tref)/2;
//
//            //- enthalpy at T
//            hb1 = species[idH2O].h(T1)*0.276+species[idN2].h(T1)*0.724;
//            hb2 = species[idH2O].h(T2)*0.276+species[idN2].h(T2)*0.724;
//            hb3 = species[idH2O].h(T3)*0.276+species[idN2].h(T3)*0.724;
//
//            //- cp at arithmetic mean temperatures
//            cpm1 = species[idH2O].cp(Tm1)*0.276+species[idN2].cp(T1)*0.724;
//            cpm2 = species[idH2O].cp(Tm2)*0.276+species[idN2].cp(T2)*0.724;
//            cpm3 = species[idH2O].cp(Tm3)*0.276+species[idN2].cp(T3)*0.724;
//
//            //- normalize
//            //  if == 0 then enthalpy of unburned equal to burned
//            normalized1 = fabs((hb1 - cpm1 * T1 - href)/href);
//            normalized2 = fabs((hb2 - cpm2 * T2 - href)/href);
//            normalized3 = fabs((hb3 - cpm3 * T3 - href)/href);
//
//            if (normalized1 < 0 && normalized3 > 0)
//            {
//                //- in interval 1 to mid point
//                T2 = T3;
//            }
//
//            if (normalized2 > 0 && normalized3 < 0)
//            {
//                //- in interval mid point to 2
//                T1 = T3;
//            }
//
//            std::cout << normalized1 << "\t" << normalized3 << "\t" << normalized2 << "\n";
//            std::cout << T1 << "\t" << T3 << "\t" << T2 << "\n";
//
//            if (fabs(normalized3) < 1e-6 || i == 30)
//            {
//                break;
//            }
//            i++; T+=100;
//        }
//        while (true);
//        std::cout<< hb3;
//        return (T3+Tref);
//    }*/
//
//    //- discret points of mixture fraction points Z [-]
//    //  TODO --- move to read afcDict function
//    scalarField discretZ
//    (
//        const normalString& fileAFCDict,
//        const scalar& Zst
//    )
//    {
//        stringField fileContent = openFile(fileAFCDict);
//        scalarField Z;
//
//        //- first entry
//        Z.push_back(0);
//
//        forAll(fileContent, line)
//        {
//
//            normalString lineContent = fileContent[line];
//
//            //- split line into array
//            std::istringstream tmp(lineContent);
//
//            stringField lineContent_
//            {
//                std::istream_iterator<std::string>{tmp},
//                std::istream_iterator<std::string>{}
//            };
//
//            if (lineContent_.size())
//            {
//                if (lineContent_[0] == "mixtureFractionPoints")
//                {
//                    scalar deltaZ = 1/stod(lineContent_[1]);
//                    scalar tmp{0};
//
//                    for(int i=1; i<atoi(lineContent_[1].c_str()); i++)
//                    {
//                        tmp += deltaZ;
//                        if (Z[i-1] < (1-Zst) && Z[i-1]+deltaZ > (1-Zst))
//                        {
//                            Z.push_back(1-Zst);
//                        }
//                        Z.push_back(tmp);
//
//                    }
//                }
//            }
//        }
//        //-last entry
//        Z.push_back(1);
//
//        return Z;
//    }
//
//    //- discret points of scalar dissipation rate [Hz]
//    //  TODO --- move to read afcDict function
//    scalarField discretChi
//    (
//        const normalString& fileAFCDict
//    )
//    {
//        stringField fileContent = openFile(fileAFCDict);
//        scalarField chi;
//        bool sDR{false};
//
//        //- first entry
//        chi.push_back(1e-8);
//
//        forAll(fileContent, line)
//        {
//
//            normalString lineContent = fileContent[line];
//
//            if (sDR)
//            {
//                removeSpace(lineContent);
//                //- dictionary end
//                if (lineContent == "}")
//                {
//                    sDR = false;
//                    break;
//                }
//
//                //- split line into array
//                //  reinitialize due to removed spaces
//                std::istringstream tmp(fileContent[line]);
//
//                stringField lineContent_
//                {
//                    std::istream_iterator<std::string>{tmp},
//                    std::istream_iterator<std::string>{}
//                };
//
//                forAll(lineContent_, num)
//                {
//                    chi.push_back(stod(lineContent_[num]));
//                }
//            }
//
//            //- remove whitespaces
//            removeSpace(lineContent);
//
//            if (lineContent == "scalarDissipationRates")
//            {
//                normalString tmp = fileContent[line+1];
//                removeSpace(tmp);
//
//                if (tmp == "{")
//                {
//                    sDR = true;
//                    line++;
//                }
//            }
//        }
//        return chi;
//    }
//
//    //- initialize mass fraction Y of fuel and oxidizer
//    void initializeY
//    (
//        const scalarField& Z_dP,
//        std::vector<std::vector<Species> >& Z
//    )
//    {
//        forAll(Z, point)
//        {
//            //- for each point loop over all species to set up the species
//            //  calculation is linear for the beginning due to no
//            //  chemical reaction. That means:
//            //  Z = 0   --> Yi = Yoi
//            //  Z = 1   --> Yi = Yfi
//            //      Yi(Z) = (Yfi - Yoi)*Z + Yoi
//            forAll(Z[point], species)
//            {
//                //- set linear mass fraction
//                Z[point][species].setY
//                (
//                    (Z[point][species].Yf() - Z[point][species].Yo()) *
//                    Z_dP[point] +
//                    Z[point][species].Yo()
//                );
//
//                //- set linear mol fraction
//                Z[point][species].setX
//                (
//                    (Z[point][species].Xf() - Z[point][species].Xo()) *
//                    Z_dP[point] +
//                    Z[point][species].Xo()
//                );
//            }
//        }
//    }
//
//
//    //- remove whitespace from string
//    void removeSpace(normalString& str)
//    {
//        str.erase
//        (
//            std::remove_if
//            (
//                str.begin(),
//                str.end(),
//                [](char c)
//                {
//                    return (c =='\r' || c =='\t' || c == ' ' || c == '\n');
//                }
//            ),
//            str.end()
//        );
//    }
//
//    //- calculate molecular weight
//    scalar calcMolecularWeight(const normalString& species)
//    {
//        //- class of molecular weights
//        Elements element;
//
//        //- molecular weight [mol/kg]
//        scalar mW{0};
//
//        //- there should be a better way to perform that stuff | simple c style
//        //  get elements and the number
//
//            normalString tmp = species;
//            normalString lP, aP, nP;
//
//            //- inert gas
//            if (tmp == "HE" || tmp == "AR")
//            {
//                mW = element.atomicWeight(tmp);
//            }
//            else
//            {
//                for (unsigned int i=0; i<tmp.length(); i++)
//                {
//                    //- actual position
//                    aP = tmp[i];
//
//                    //- next position
//                    if (i+1 <= tmp.length())
//                    {
//                        nP = tmp[i+1];
//                    }
//                    else
//                    {
//                        nP = aP;
//                    }
//
//                    //- if aP is no number
//                    //- TODO add two numbered molecules like C4H12 <-- 12
//                    if (aP.find_first_of("0123456789") == std::string::npos)
//                    {
//                        if
//                        (
//                            nP.find_first_of("0123456789") != std::string::npos
//                        )
//                        {
//                            mW += element.atomicWeight(aP)
//                                * stod(nP);
//
//                            i++;
//                        }
//                        //- only one element
//                        else
//                        {
//                            mW += element.atomicWeight(aP);
//                        }
//                    }
//                }
//            }
//
//            //- return the molecular weight of species
//            return mW;
//    }
//
//    //- summary
//    void summary
//    (
//        const scalar& hf_a,
//        const scalar& ho_a,
//        const scalar& Zst,
//        const scalarField& chi_dP,
//        //const scalar& Tst_a,
//        const std::vector<Species>& species
//    )
//    {
//        scalar Tf{0};
//        scalar To{0};
//
//        forAll(species, id)
//        {
//            if (species[id].fuel())
//            {
//                Tf = species[id].Tf();
//            }
//            if (species[id].oxidizer())
//            {
//                To = species[id].To();
//            }
//        }
//
//        std::cout<< std::setiosflags(std::ios::left)
//                 << "\n\nSummary of important data\n"
//                 <<"------------------------------------------------------\n"
//                 << std::setw(25)<< "Temperature fuel: " << Tf << "\t [K]\n"
//                 << std::setw(25)<< "Enthalpy fuel: " << hf_a << "\t [J/kg]\n"
//                 << std::setw(25)<< "Temperature oxidizer: " << To << "\t [K]\n"
//                 << std::setw(25)<< "Enthalpy oxidizer: " << ho_a << "\t [J/kg]\n"
//                 << std::setw(25)<< "Stochiometric: " << Zst << "\t [-]\n"
//                 << "Discrete scalar dissipation rates: " << chi_dP.size() << "\t [-]\n"
//                 << "Range of scalar dissipation rates: " << chi_dP[0] << " - " << chi_dP[chi_dP.size()-1] << "\t [Hz]\n"
////                 << "Adiabatic flame temperature: " << Tst_a << "[K]\n"
//                 << "\n\n";
//    }
