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
#ifndef transportProperties_hpp
#define transportProperties_hpp


class TransportProperties
{
    public:

        //- constructor
        TransportProperties();

        //- destructor
        ~TransportProperties();


    public:

    //- functions
    //-----------

        //- set geometry;
        void setGeometry( const int& );

        //- get geometry
        int getGeometry() const;

        //- set Lennard-Jones potential
        void setLJPotential( const float& );

        //- get Lennard-Jones potential
        float getLJPotential() const;

        //- set Lennard-Jones collision diameter
        void setLJCollisionDiameter( const float& );

        //- get Lennard-Jones collision diameter
        float getLJCollisionDiameter() const;

        //- set dipole moment in Debye
        void setDipoleMoment( const float& );

        //- get dipole moment in Debye
        float getDipoleMoment() const;

        //- set polarizability in Angstroms
        void setPolarizability( const float& );

        //- get polarizability in Angstroms
        float getPolarizability() const;

        //- set rotational relaxation collision number
        void setRotRelaxCollNo( const float& );

        //- get rotational relaxation collision number
        float getRotRelaxCollNo() const;

        //- get values of species
        void showValues() const;


    private:

    //- variables
    //-----------

        //- geometry
        int geometry{0};

        //- the Lennard-Jones potential well depth ε/k_B in [K]
        float eps_kB{0.};

        //- the Lennard-Jones collision diameter σ in Angstroms == [pm] (pico meter)
        float sigma{0.};

        //- the dipole moment μ in Debye. Note: a Debye is 3.33564*10^-10 [Cm]
        float mu{0.};

        //- the polarizability α in cubic Angstroms == [pm]
        float alpha{0.};

        //- the rotational relaxation collision number Zrot at 298K
        float Zrot{0.};
};


#endif // transportProperties_hpp
