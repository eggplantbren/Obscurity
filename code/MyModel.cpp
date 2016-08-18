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
,image(ni, std::vector<double>(nj, 0.0))
,x(nj)
,y(ni)
{
    if(std::abs((dx - dy)/dx) > 1E-6)
        throw std::domain_error("Non-square pixels.");

    for(size_t i=0; i<x.size(); ++i)
        x[i] = x_min + (i + 0.5)*dx;

    for(size_t i=0; i<y.size(); ++i)
        y[i] = y_max - (i + 0.5)*dy;
}

void MyModel::from_prior(DNest4::RNG& rng)
{
    blobs.from_prior(rng);
    x0 = x_min - (x_max - x_min) + (x_max - x_min)*rng.rand();

    DNest4::Cauchy c(0.0, data.get_t_range());
    timescale = std::abs(c.generate(rng));

    calculate_image();
}

void MyModel::calculate_image()
{
    double star_density = 1.0/M_PI;
    double rsq;
    double f;
    for(size_t i=0; i<image.size(); ++i)
    {
        for(size_t j=0; j<image.size(); ++j)
        {
            rsq = x[j]*x[j] + y[i]*y[i];
            if(rsq < 1.0)
            {
                f = blobs.evaluate(x[j], y[i]);
                if(f < 10.0)
                    image[i][j] = star_density*exp(-f);
                else
                    image[i][j] = 0.0;
            }
            else
                image[i][j] = 0.0;
        }
    }
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

    calculate_image();

	return logH;
}

double MyModel::log_likelihood() const
{
	double logL = 0.0;
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

