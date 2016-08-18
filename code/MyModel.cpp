#include "MyModel.h"
#include "DNest4/code/DNest4.h"
#include <stdexcept>

namespace Obscurity
{

// STATIC STUFF
Data MyModel::data;

const double MyModel::dx = (MyModel::x_max - MyModel::x_min)/MyModel::nj;
const double MyModel::dy = (MyModel::y_max - MyModel::y_min)/MyModel::ni;

void MyModel::load_data(const char* filename)
{
    data.load(filename);
}

// CONSTRUCTOR AND MEMBER FUNCTIONS
MyModel::MyModel()
:blobs(x_min, x_max, y_min, y_max)
{
    if(std::abs((dx - dy)/dx) > 1E-6)
        throw std::domain_error("Non-square pixels.");
}

void MyModel::from_prior(DNest4::RNG& rng)
{
    blobs.from_prior(rng);
    x0 = x_min - (x_max - x_min) + (x_max - x_min)*rng.rand();

    DNest4::Cauchy c(0.0, data.get_t_range());
    timescale = std::abs(c.generate(rng));
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
        x0 += (x_max - x_min)*rng.randh();
        DNest4::wrap(x0, x_min - (x_max - x_min), 0.0);
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

    const auto& t = data.get_t();
    const auto& y = data.get_y();
    const auto& sig = data.get_sig();

	return logL;
}

void MyModel::print(std::ostream& out) const
{
    out<<x0<<' '<<timescale<<' ';

//    for(size_t i=0; i<image.size(); ++i)
//        for(size_t j=0; j<image.size(); ++j)
//            out<<image[i][j]<<' ';
}

std::string MyModel::description() const
{
	return std::string("");
}

} // namespace Obscurity

