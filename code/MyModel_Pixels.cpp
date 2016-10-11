#include "MyModel_Pixels.h"
#include "DNest4/code/DNest4.h"
#include <stdexcept>
#include <armadillo>

namespace Obscurity
{

// CONSTRUCTOR AND MEMBER FUNCTIONS
MyModel_Pixels::MyModel_Pixels()
:n(ni, nj)
,obscurer_map(ni, nj)
,convolved(ni, nj)
{

}

void MyModel_Pixels::from_prior(DNest4::RNG& rng)
{
    for(size_t j=0; j<nj; ++j)
        for(size_t i=0; i<ni; ++i)
            n(i, j) = rng.randn();

    DNest4::Cauchy c1(0.0, 1.0);
    x0 = -std::abs(c1.generate(rng));

    DNest4::Cauchy c2(0.0, 0.1*data.get_t_range());
    timescale = std::abs(c2.generate(rng));

    calculate_obscurer_map();
}

double MyModel_Pixels::calculate_total_flux(double time) const
{
    double speed = 1.0/timescale;
    double offset = x0 + speed*time;

    int j = (int)floor((offset - (x_min + 0.5*dx))/dx);
    return convolved(ni/2, j%nj);
}

double MyModel_Pixels::perturb(DNest4::RNG& rng)
{
    double logH = 0.0;

    if(rng.rand() <= 0.7)
    {

    }
    else if(rng.rand() <= 0.5)
    {
        DNest4::Cauchy c(0.0, 1.0);
        logH += c.perturb(x0, rng);
        x0 = -std::abs(x0);
    }
    else
    {
        DNest4::Cauchy c(0.0, 0.1*data.get_t_range());
        logH += c.perturb(timescale, rng);
        timescale = std::abs(timescale);
    }
    return logH;
}

double MyModel_Pixels::log_likelihood() const
{
    double logL = 0.0;

    const auto& t = data.get_t();
    const auto& y = data.get_y();
    const auto& sig = data.get_sig();

    double model_prediction;
    for(size_t i=0; i<t.size(); ++i)
    {
        model_prediction = calculate_total_flux(t[i]);
        logL += -0.5*log(2*M_PI) - log(sig[i])
                    - 0.5*pow((y[i] - model_prediction)/sig[i], 2);
    }

    return logL;
}

void MyModel_Pixels::calculate_obscurer_map()
{
    // Obscurer image
    for(size_t j=0; j<nj; ++j)
        for(size_t i=0; i<ni; ++i)
            obscurer_map(i, j) = 0.0;

    for(size_t j=0; j<nj; ++j)
        for(size_t i=0; i<ni; ++i)
            obscurer_map(i, j) = exp(-obscurer_map(i, j));

    // FFT of obscurer_map
    arma::cx_mat A = arma::fft2(obscurer_map);

    // (FFT of obscurer map)*(FFT of star map)
    for(size_t j=0; j<nj; ++j)
        for(size_t i=0; i<ni; ++i)
            A(i, j) *= fft_of_star(i, j);

    // obscurer map convolved with star map
    A = arma::ifft2(A);
    for(size_t i=0; i<ni; ++i)
        for(size_t j=0; j<nj; ++j)
            convolved(i, j) = A(i, j).real();
}

void MyModel_Pixels::print(std::ostream& out) const
{
    const auto& t = data.get_t();
    for(size_t i=0; i<t.size(); ++i)
        out<<calculate_total_flux(t[i])<<' ';

    for(size_t i=0; i<ni; ++i)
        for(size_t j=0; j<nj; ++j)
            out<<star(i, j)<<' ';

    for(size_t i=0; i<ni; ++i)
        for(size_t j=0; j<nj; ++j)
            out<<obscurer_map(i, j)<<' ';
}

std::string MyModel_Pixels::description() const
{
    return std::string("");
}

/* STATIC STUFF */

Data                  MyModel_Pixels::data;
arma::mat             MyModel_Pixels::star(MyModel_Pixels::ni, MyModel_Pixels::nj);
arma::cx_mat          MyModel_Pixels::fft_of_star(MyModel_Pixels::ni, MyModel_Pixels::nj);
std::vector<double>   MyModel_Pixels::x(MyModel_Pixels::nj);
std::vector<double>   MyModel_Pixels::y(MyModel_Pixels::ni);

void MyModel_Pixels::initialise()
{
    for(size_t i=0; i<ni; ++i)
        y[i] = y_max - (i + 0.5)*dy;

    for(size_t j=0; j<nj; ++j)
        x[j] = x_min + (j + 0.5)*dx;

    double rsq;

    double limb_darkening_coefficient = 1.0;

    // Argh column-major order
    // Star image
    double tot = 0.0;
    for(size_t j=0; j<nj; ++j)
    {
        for(size_t i=0; i<ni; ++i)
        {
            rsq = x[j]*x[j] + y[i]*y[i];
            if(rsq < 1.0)
            {
                star(i, j) = 1.0 -
                    limb_darkening_coefficient*(1.0 - sqrt(1.0 - rsq));
            }
            else
                star(i, j) = 0.0;
            tot += star(i, j);
        }
    }
    for(size_t j=0; j<nj; ++j)
        for(size_t i=0; i<ni; ++i)
            star(i, j) /= tot;

    arma::mat star2 = star;
    int m, n;
    for(int i=0; i<(int)ni; i++)
    {
        m = DNest4::mod(i - (int)ni/2, (int)ni);
        for(int j=0; j<(int)nj; j++)
        {
            n = DNest4::mod(j - (int)nj/2, (int)nj);
            star2(m, n) = star(i, j);
        }
    }

    fft_of_star = arma::fft2(star2);
}

void MyModel_Pixels::load_data(const char* filename)
{
    data.load(filename);
}

} // namespace Obscurity

