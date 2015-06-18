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
#include "chemistry.hpp"


Chemistry::Chemistry() : r_(0)
{
}


Chemistry::~Chemistry()
{
}


void Chemistry::set_kf(const bool a)
{
    //- check the matrix
    //  if matrix elements equal or higher than
    //  the number reaction, there is sth. wron
    if (fb_.size() >= r_)
    {
        std::cerr<< "\n ++ ERROR: fb matrix rows >= reactions "
                 << fb_.size() << " >= " << r_
                 << "\n ++ Error occur in file " << __FILE__
                 << " line " << __LINE__ << std::endl;
        std::terminate();
    }

    fb_.push_back(std::vector<double>(0));

    //- add two new columns
    fb_[r_-1].push_back(0);
    fb_[r_-1].push_back(0);

    //- forward reaction used
    if (a)
    {
        fb_[r_-1][0] = 1;
    }
    //- forward reaction NOT used
    else
    {
        fb_[r_-1][0] = 0;
    }
}

void Chemistry::set_kb(const bool a)
{
    //- check the matrix
    //  Due to the fact that set_kb is after set_kf
    //  the number reactions has to be equal to the
    //  matrix size
    if (fb_.size() != r_)
    {
        std::cerr<< "\n ++ ERROR: fb matrix rows != reactions "
                 << fb_.size() << " != " << r_
                 << "\n ++ Error occur in file " << __FILE__
                 << " line " << __LINE__ << std::endl;
        std::terminate();
    }

    //- backward reaction used
    if (a)
    {
        fb_[r_-1][1] = 1;
    }
    //- backward reaction NOT used
    else
    {
        fb_[r_-1][1] = 0;
    }
}


unsigned int Chemistry::fb
(
    const unsigned int r,
    const unsigned int i
) const
{
    //- return (f)orward and (b)ackward coeffs
    //  i = 0   (f)
    //  i = 1   (b)
    return fb_[r][i];
}


void Chemistry::setElementarReaction(const normalString& formula)
{
    //- check the vector
    //  if matrix elements equal or higher than
    //  the number reaction, there is sth. wrong
    if (elementarReaction_.size() >= r_)
    {
        std::cerr<< "\n ++ ERROR: elementarReaction vector > reactions... "
                 << "\n ++ Error occur in file " << __FILE__
                 << " line " << __LINE__ << std::endl;
        std::terminate();
    }

    elementarReaction_.push_back(formula);
}


normalString Chemistry::elementarReaction(const unsigned int& i) const
{
    return elementarReaction_[i];
}


void Chemistry::increment_r()
{
    r_++;
}


void Chemistry::decrement_r()
{
    r_--;
}


int Chemistry::r() const
{
    return r_;
}


void Chemistry::setTROE(const bool a)
{
    //- no check of vector size because
    //  if there should be a problem we will get the
    //  trouble in the kf_ kb_ stuff
    if (TROE_.size() < r_)
    {
        TROE_.push_back(0);
    }

    //- TROE used
    if (a)
    {
        TROE_[r_-1] = 1;
    }
    //- TROE NOT used
    else
    {
        TROE_[r_-1] = 0;
    }
}


void Chemistry::setLOW(const bool a)
{
    //- no check of vector size because
    //  if there should be a problem we will get the
    //  trouble in the kf_ kb_ stuff
    if (LOW_.size() < r_)
    {
        LOW_.push_back(0);
    }

    //- TROE used
    if (a)
    {
        LOW_[r_-1] = 1;
    }
    //- TROE NOT used
    else
    {
        LOW_[r_-1] = 0;
    }
}


unsigned int Chemistry::TROE(const unsigned int r) const
{
    return TROE_[r];
}


unsigned int Chemistry::LOW(const unsigned int r) const
{
    return LOW_[r];
}


void Chemistry::setArrheniusCoeffs
(
    const scalar& A0,
    const scalar& n,
    const scalar& E0
)
{
    //- no check of vector size because
    //  if there should be a problem we will get the
    //  trouble in the kf_ kb_ stuff
    std::cout << arrheniusCoeffs_.size() << " - " << r_ << std::endl;
    if (arrheniusCoeffs_.size() < r_)
    {
        arrheniusCoeffs_.push_back(std::vector<double>(0));
    }

    arrheniusCoeffs_.push_back(std::vector<double>(0));

    //- add three new columns
    arrheniusCoeffs_[r_-1].push_back(A0);
    arrheniusCoeffs_[r_-1].push_back(n);
    arrheniusCoeffs_[r_-1].push_back(E0);
}


scalar Chemistry::arrheniusCoeffs
(
    const unsigned int& r,
    const unsigned int i
) const
{
    //- return coeffs
    //  i = 0 -> A0
    //  i = 1 -> n
    //  i = 2 -> Ea
    return arrheniusCoeffs_[r][i];
}
