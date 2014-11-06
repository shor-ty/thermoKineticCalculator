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
#include <iostream>
//- user def. headers
#include "transportProperties.hpp"

//- constructor
TransportProperties::TransportProperties()
{
}

//- destructor
TransportProperties::~TransportProperties()
{
}

//- functions
//-----------

    //- set geometry factor
    void TransportProperties::setGeometry( const int& speciesGeometry )
    {
        geometry = speciesGeometry;
    }

    //- set potential
    void TransportProperties::setLJPotential( const float& speciesLJPotential )
    {
        eps_kB = speciesLJPotential;
    }

    //- set collision diameter
    void TransportProperties::setLJCollisionDiameter( const float& diameter )
    {
        sigma = diameter;
    }

    //- set dipole momentum
    void TransportProperties::setDipoleMoment( const float& dipoleMoment )
    {
        mu = dipoleMoment;
    }

    //- set polarizability
    void TransportProperties::setPolarizability( const float& polarizability )
    {
        alpha = polarizability;
    }

    //- set rotation relaxation collision number
    void TransportProperties::setRotRelaxCollNo( const float& rotRelaxCollNo )
    {
        Zrot = rotRelaxCollNo;
    }

    //- show values
    void TransportProperties::showValues() const
    {
        std::cout << " - Transport properties:\n\n";
        std::cout << "   + Geometry factor:                        " << geometry << "\n";
        std::cout << "   + Dipole moment in Debye:                 " << mu << "\n";
        std::cout << "   + Lennard-Jones potential:                " << eps_kB << "\n";
        std::cout << "   + Polarizability in Angstroms:            " << alpha << "\n";
        std::cout << "   + Lennard-Jones collision diameter:       " << sigma<< "\n";
        std::cout << "   + Rotational relaxation collision number: " << alpha << std::endl;
    }


//- functions
//-----------

    //- get geometry
    int TransportProperties::getGeometry() const
    {
        return geometry;
    }

    //- get potential
    float TransportProperties::getLJPotential() const
    {
        return eps_kB;
    }

    //- get collision diameter
    float TransportProperties::getLJCollisionDiameter() const
    {
        return sigma;
    }

    //- get dipole momentum
    float TransportProperties::getDipoleMoment() const
    {
        return mu;
    }

    //- get polarizability diameter
    float TransportProperties::getPolarizability() const
    {
        return alpha;
    }

    //- get rotation relaxation collision no
    float TransportProperties::getRotRelaxCollNo() const
    {
        return Zrot;
    }
