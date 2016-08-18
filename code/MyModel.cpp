#include "MyModel.h"
#include "DNest4/code/DNest4.h"
#include <stdexcept>
#include <armadillo>

namespace Obscurity
{

// CONSTRUCTOR AND MEMBER FUNCTIONS
MyModel::MyModel()
:blobs(4, 100, false, MyConditionalPrior(), DNest4::PriorType::log_uniform)
{

}

void MyModel::from_prior(DNest4::RNG& rng)
{
    blobs.from_prior(rng);

    DNest4::Cauchy c1(0.0, 1.0);
    x0 = -std::abs(c1.generate(rng));

    DNest4::Cauchy c2(0.0, data.get_t_range());
    timescale = std::abs(c2.generate(rng));
}

double MyModel::calculate_total_flux(double time) const
{
    double speed = 1.0/timescale;
    double offset = x0 + speed*time;

    return 0.0;
}

double MyModel::perturb(DNest4::RNG& rng)
{
	double logH = 0.0;

    if(rng.rand() <= 0.7)
    {
        logH += blobs.perturb(rng);
    }
    else if(rng.rand() <= 0.5)
    {
        DNest4::Cauchy c(0.0, 1.0);
        logH += c.perturb(x0, rng);
        x0 = -std::abs(x0);
    }
    else
    {
        DNest4::Cauchy c(0.0, data.get_t_range());
        logH += c.perturb(timescale, rng);
        timescale = std::abs(timescale);
    }
	return logH;
}

double MyModel::log_likelihood() const
{
	double logL = 0.0;

//    const auto& t = data.get_t();
//    const auto& y = data.get_y();
//    const auto& sig = data.get_sig();

    arma::mat star(ni, nj);
    arma::mat obscurer(ni, nj);



    // Obscurer image
    for(size_t j=0; j<nj; ++j)
        for(size_t i=0; i<ni; ++i)
            obscurer(i, j) = 0.0;

    const auto& blobs_params = blobs.get_components();
    int i_min, i_max, j_min, j_max;
    for(const auto& blob_params: blobs_params)
    {
        // Determine square patch of image to loop over
        i_min = (int)floor(((y_max - 0.5*dy) - (blob_params[1] + blob_params[3]))/dy);
        i_max = (int)floor(((y_max - 0.5*dy) - (blob_params[1] - blob_params[3]))/dy);
        j_min = (int)floor(((blob_params[0] - blob_params[3]) - (x_min + 0.5*dx))/dx);
        j_max = (int)floor(((blob_params[0] + blob_params[3]) - (x_min + 0.5*dx))/dx);

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
                obscurer(i, j) += evaluate_blob(blob_params, x[j], y[i]);
    }

	return logL;
}

void MyModel::print(std::ostream& out) const
{
//    if(std::abs(dx - dy) > 1E-6*sqrt(dx*dy))
//        throw std::domain_error("Non-square pixels.");

}

std::string MyModel::description() const
{
	return std::string("");
}

/* STATIC STUFF */

Data                  MyModel::data;
arma::mat             MyModel::star(MyModel::ni, MyModel::nj);
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
    for(size_t j=0; j<nj; ++j)
    {
        for(size_t i=0; i<ni; ++i)
        {
            rsq = x[j]*x[j] + y[i]*y[i];
            if(rsq <= 1.0)
                star(i, j) = 1.0;
            else
                star(i, j) = 0.0;
        }
    }
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
    double widthsq = blob_params[3]*blob_params[3];

    if(rsq >= widthsq)
        return 0.0;

    return C*blob_params[2]*(1.0 - rsq/widthsq)/widthsq;
}


} // namespace Obscurity

