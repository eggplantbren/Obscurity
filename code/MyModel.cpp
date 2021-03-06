#include "MyModel.h"
#include "DNest4/code/DNest4.h"
#include <stdexcept>
#include <cmath>
#include <armadillo>
#include <sstream>

namespace Obscurity
{

// CONSTRUCTOR AND MEMBER FUNCTIONS
MyModel::MyModel()
:blobs(4, 100, false, MyConditionalPrior(), DNest4::PriorType::log_uniform)
,obscurer_map(ni, nj)
,convolved(ni, nj)
,star(ni, nj)
,fft_of_star(ni, nj)
{

}

void MyModel::from_prior(DNest4::RNG& rng)
{
    blobs.from_prior(rng);

    DNest4::Cauchy c1(0.0, 1.0);
    x0 = -std::abs(c1.generate(rng));

    DNest4::Cauchy c2(0.0, 0.1*data.get_t_range());
    timescale = std::abs(c2.generate(rng));

    // Limb darkening prior
    // p(x) \propto 0.7752 + 1.4379*x^2.3369
    auto p = [](double x) { return 0.7752 + 1.4379*pow(x, 2.3369); };
    double p_max = p(1.0);
    while(true)
    {
        limb_darkening = rng.rand(); 
        if(rng.rand() <= p(limb_darkening)/p_max)
            break;
    }

    DNest4::Cauchy c3(0.0, 5.0);
    magnitude = c3.generate(rng);

    noise_boost = rng.rand();
    nu = exp(log(1.0) + log(1E3)*rng.rand());

    calculate_star();
    calculate_obscurer_map();
}

double MyModel::calculate_total_flux(double time) const
{
    double speed = 1.0/timescale;
    double offset = x0 + speed*time;

    int j = (int)floor((offset - (x_min + 0.5*dx))/dx);
    return convolved(ni/2, j%nj);
}

double MyModel::perturb(DNest4::RNG& rng)
{
    double logH = 0.0;

    if(rng.rand() <= 0.7)
    {
        logH += blobs.perturb(rng);
        if(blobs.components_changed())
            calculate_obscurer_map();
    }
    else
    {
        int which = rng.rand_int(6);
        if(which == 0)
        {
            DNest4::Cauchy c(0.0, 1.0);
            logH += c.perturb(x0, rng);
            x0 = -std::abs(x0);
        }
        else if(which == 1)
        {
            DNest4::Cauchy c(0.0, 0.1*data.get_t_range());
            logH += c.perturb(timescale, rng);
            timescale = std::abs(timescale);
        }
        else if(which == 2)
        {
            // Limb darkening prior density
            auto p = [](double x) { return 0.7752 + 1.4379*pow(x, 2.3369); };
            logH -= log(p(limb_darkening));
            limb_darkening += rng.randh();
            DNest4::wrap(limb_darkening, 0.0, 1.0);
            logH += log(p(limb_darkening));

            // Recalculate star
            calculate_star();
        }
        else if(which == 3)
        {
            DNest4::Cauchy c3(0.0, 5.0);
            logH += c3.perturb(magnitude, rng);
        }
        else if(which == 4)
        {
            noise_boost += rng.randh();
            DNest4::wrap(noise_boost, 0.0, 1.0);
        }
        else
        {
            nu = log(nu);
            nu += log(1E3)*rng.randh();
            DNest4::wrap(nu, log(1.0), log(1E3));
            nu = exp(nu);
        }
    }
    return logH;
}

double MyModel::log_likelihood() const
{
    double logL = 0.0;

    const auto& t = data.get_t();
    const auto& y = data.get_y();
    const auto& sig = data.get_sig();

    // Transform from latent U(0, 1) to actual noise-boost factor
    double K = (noise_boost < 0.5)?
               (1.0):
               (exp(log(1E3)*2*(noise_boost - 0.5)));

    // Compute intensive constant factors outside the loop
    double C = std::lgamma(0.5*(nu + 1.0)) - std::lgamma(0.5*nu)
                - 0.5*log(nu*M_PI);

    double model_prediction, scale;
    for(size_t i=0; i<t.size(); ++i)
    {
        model_prediction = magnitude - 2.5*log10(calculate_total_flux(t[i]));
        scale = K*sig[i];
        logL += C - log(scale) - 0.5*(nu+1)*
                        (1.0 + pow(y[i] - model_prediction, 2)/scale/scale/nu);
    }

    return logL;
}

void MyModel::calculate_obscurer_map()
{
    // Obscurer image
    for(size_t j=0; j<nj; ++j)
        for(size_t i=0; i<ni; ++i)
            obscurer_map(i, j) = 0.0;

    const auto& blobs_params = blobs.get_components();
    int i_min, i_max, j_min, j_max;
    double width;

    // Center of mass of blobs
    auto com = com_blobs();

    for(const auto& blob_params: blobs_params)
    {
        // Use position relative to center of mass
        std::vector<double> centered = blob_params;
        centered[0] -= std::get<0>(com);
        centered[1] -= std::get<1>(com);    

        width = centered[3]*LL;

        // Determine square patch of image to loop over
        i_min = (int)floor(((y_max - 0.5*dy) - (centered[1] + width))/dy);
        i_max = (int)floor(((y_max - 0.5*dy) - (centered[1] - width))/dy);
        j_min = (int)floor(((centered[0] - width) - (x_min + 0.5*dx))/dx);
        j_max = (int)floor(((centered[0] + width) - (x_min + 0.5*dx))/dx);

        if(i_min < 0)
            i_min = 0;
        if(i_max < 0)
            i_max = 0;
        if(j_min < 0)
            j_min = 0;
        if(j_max < 0)
            j_max = 0;

        if(i_min >= (int)ni)
            i_min = ni - 1;
        if(i_max >= (int)ni)
            i_max = ni - 1;
        if(j_min >= (int)nj)
            j_min = nj - 1;
        if(j_max >= (int)nj)
            j_max = nj - 1;

        for(int j=j_min; j<=j_max; ++j)
            for(int i=i_min; i<=i_max; ++i)
                obscurer_map(i, j) += evaluate_blob(centered, x[j], y[i]);
    }
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

void MyModel::calculate_star()
{
    double rsq;

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
                    limb_darkening*(1.0 - sqrt(1.0 - rsq));
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

void MyModel::print(std::ostream& out) const
{
    // Transform from latent U(0, 1) to actual noise-boost factor
    double K = (noise_boost < 0.5)?
               (1.0):
               (exp(log(1E3)*2*(noise_boost - 0.5)));

    out<<x0<<' '<<timescale<<' '<<limb_darkening<<' ';
    out<<K<<' '<<nu<<' ';
    blobs.print(out);

    const auto& t = data.get_t();
    for(size_t i=0; i<t.size(); ++i)
        out<<(magnitude - 2.5*log10(calculate_total_flux(t[i])))<<' ';

    for(size_t i=0; i<ni; ++i)
        for(size_t j=0; j<nj; ++j)
            out<<star(i, j)<<' ';

    for(size_t i=0; i<ni; ++i)
        for(size_t j=0; j<nj; ++j)
            out<<obscurer_map(i, j)<<' ';
}

std::string MyModel::description() const
{
    std::stringstream s;
    s<<"x0, timescale, limb_darkening, K, nu, ";

    s<<"dim_blobs, max_num_blobs, sigma, mu_mass, mu_width, ";
    s<<"num_blobs, ";
    for(int i=0; i<100; ++i)
        s<<"xc["<<i<<"], ";
    for(int i=0; i<100; ++i)
        s<<"yc["<<i<<"], ";
    for(int i=0; i<100; ++i)
        s<<"mass["<<i<<"], ";
    for(int i=0; i<100; ++i)
        s<<"width["<<i<<"], ";

    const auto& t = data.get_t();
    for(size_t i=0; i<t.size(); ++i)
        s<<"model_prediction["<<i<<"], ";

    for(size_t i=0; i<ni; ++i)
        for(size_t j=0; j<nj; ++j)
            s<<"star("<<i<<", "<<j<<"), ";

    for(size_t i=0; i<ni; ++i)
        for(size_t j=0; j<nj; ++j)
            s<<"obscurer_map("<<i<<", "<<j<<"), ";

    return s.str();
}

/* STATIC STUFF */

Data                  MyModel::data;
std::vector<double>   MyModel::x(MyModel::nj);
std::vector<double>   MyModel::y(MyModel::ni);

void MyModel::initialise()
{
    for(size_t i=0; i<ni; ++i)
        y[i] = y_max - (i + 0.5)*dy;

    for(size_t j=0; j<nj; ++j)
        x[j] = x_min + (j + 0.5)*dx;
}

void MyModel::load_data(const char* filename)
{
    data.load(filename);
}

std::tuple<double, double> MyModel::com_blobs() const
{
    double com_x = 0.0;
    double com_y = 0.0;
    double mtot = 0.0;
    const auto& blobs_params = blobs.get_components();

    double xc, yc, mass;
    for(size_t i=0; i<blobs_params.size(); ++i)
    {
        xc = blobs_params[i][0];
        yc = blobs_params[i][1];
        mass = blobs_params[i][2];

        com_x += mass*xc;
        com_y += mass*yc;
        mtot += mass;
    }

    if(blobs_params.size() != 0)
    {
        com_x /= mtot;
        com_y /= mtot;
    }

    return std::tuple<double, double>(com_x, com_y);
}

double MyModel::evaluate_blob(const std::vector<double>& blob_params,
                                                        double x, double y)
{
    static constexpr double C = 2/M_PI;

    double rsq = pow(x - blob_params[0], 2) + pow(y - blob_params[1], 2);
    double widthsq = blob_params[3]*blob_params[3]*LL*LL;

    if(rsq >= widthsq)
        return 0.0;

    return C*blob_params[2]*(1.0 - rsq/widthsq)/widthsq;
}


} // namespace Obscurity

