#include "MyModel.h"
#include "DNest4/code/DNest4.h"
#include <stdexcept>
#include <armadillo>

namespace Obscurity
{

// CONSTRUCTOR AND MEMBER FUNCTIONS
MyModel::MyModel()
:blobs(4, 100, false, MyConditionalPrior(), DNest4::PriorType::log_uniform)
,obscurer_map(ni, nj)
,convolved(ni, nj)
{

}

void MyModel::from_prior(DNest4::RNG& rng)
{
    blobs.from_prior(rng);

    DNest4::Cauchy c1(0.0, 1.0);
    x0 = -std::abs(c1.generate(rng));

    DNest4::Cauchy c2(0.0, 0.1*data.get_t_range());
    timescale = std::abs(c2.generate(rng));

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

double MyModel::log_likelihood() const
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

void MyModel::calculate_obscurer_map()
{
    // Obscurer image
    for(size_t j=0; j<nj; ++j)
        for(size_t i=0; i<ni; ++i)
            obscurer_map(i, j) = 0.0;

    const auto& blobs_params = blobs.get_components();
    int i_min, i_max, j_min, j_max;
    double width;

    for(const auto& blob_params: blobs_params)
    {
        width = blob_params[3]*LL;

        // Determine square patch of image to loop over
        i_min = (int)floor(((y_max - 0.5*dy) - (blob_params[1] + width))/dy);
        i_max = (int)floor(((y_max - 0.5*dy) - (blob_params[1] - width))/dy);
        j_min = (int)floor(((blob_params[0] - width) - (x_min + 0.5*dx))/dx);
        j_max = (int)floor(((blob_params[0] + width) - (x_min + 0.5*dx))/dx);

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
                obscurer_map(i, j) += evaluate_blob(blob_params, x[j], y[i]);
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

void MyModel::print(std::ostream& out) const
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

std::string MyModel::description() const
{
	return std::string("");
}

/* STATIC STUFF */

Data                  MyModel::data;
arma::mat             MyModel::star(MyModel::ni, MyModel::nj);
arma::cx_mat          MyModel::fft_of_star(MyModel::ni, MyModel::nj);
std::vector<double>   MyModel::x(MyModel::nj);
std::vector<double>   MyModel::y(MyModel::ni);

void MyModel::initialise()
{
    for(size_t i=0; i<ni; ++i)
        y[i] = y_max - (i + 0.5)*dy;

    for(size_t j=0; j<nj; ++j)
        x[j] = x_min + (j + 0.5)*dx;

    double rsq;

    // Argh column-major order
    // Star image
    double tot = 0.0;
    for(size_t j=0; j<nj; ++j)
    {
        for(size_t i=0; i<ni; ++i)
        {
            rsq = x[j]*x[j] + y[i]*y[i];
            if(rsq <= 1.0)
                star(i, j) = 1.0;
            else
                star(i, j) = 0.0;
            tot += star(i, j);
        }
    }
    for(size_t j=0; j<nj; ++j)
        for(size_t i=0; i<ni; ++i)
            star(i, j) /= tot;

    fft_of_star = arma::fft2(star);
}

void MyModel::load_data(const char* filename)
{
    data.load(filename);
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

