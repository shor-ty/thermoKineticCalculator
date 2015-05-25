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
#include "reactions.hpp"


Reactions::Reactions() :

    //- initialisation step
    kf_(false), kb_(false),

    LOW{false}, Ea_LOW(0), A_LOW(0), b_LOW(0),

    TROE{false}, a{0}, Ts{0}, Tss{0}, Tsss{0},

    Ea(0), A(0), b(0)
{
}


Reactions::~Reactions()
{
}


void Reactions::set_kf()
{
    kf_ = true;
}


void Reactions::set_kb()
{
    kb_ = true;
}


bool Reactions::kf() const
{
    return kf_;
}


bool Reactions::kb() const
{
    return kb_;
}


void Reactions::setElementarReaction(const normalString& elementarReaction)
{
    elementarReaction_ = elementarReaction;
}


normalString Reactions::elementarReaction() const
{
    return elementarReaction_;
}


void Reactions::setArrheniusCoeffs
(
    const scalar& A_,
    const scalar& b_,
    const scalar& Ea_
)
{
    if (LOW)
    {
        A_LOW = A_;
        b_LOW = b_;
        Ea_LOW = Ea_;
    }
    else
    {
        A = A_;
        b = b_;
        Ea = Ea_;
    }
}


void Reactions::setTROECoeffs
(
    const scalar& a_,
    const scalar& Tsss_,
    const scalar& Ts_,
    const scalar& Tss_
)
{
    a = a_;
    Ts = Ts_;
    Tss = Tss_;
    Tsss = Tsss_;
}


void Reactions::setLOW()
{
    LOW = true;
}


bool Reactions::statusLOW() const
{
    return LOW;
}


void Reactions::setTROE()
{
    TROE = true;
}


bool Reactions::statusTROE() const
{
    return TROE;
}


scalar Reactions::TROE_a() const
{
    return a;
}

scalar Reactions::TROE_Ts() const
{
    return Ts;
}

scalar Reactions::TROE_Tss() const
{
    return Tss;
}

scalar Reactions::TROE_Tsss() const
{
    return Tsss;
}



scalar Reactions::arrhenius_A() const
{
    return A;
}


scalar Reactions::arrhenius_b() const
{
    return b;
}


scalar Reactions::arrhenius_Ea() const
{
    return Ea;
}


scalar Reactions::arrhenius_A_LOW() const
{
    return A_LOW;
}


scalar Reactions::arrhenius_b_LOW() const
{
    return b_LOW;
}


scalar Reactions::arrhenius_Ea_LOW() const
{
    return Ea_LOW;
}

