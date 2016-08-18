#include "MyModel.h"
#include "DNest4/code/DNest4.h"
#include <stdexcept>

namespace Obscurity
{

// STATIC STUFF
Data MyModel::data;

void MyModel::load_data(const char* filename)
{
    data.load(filename);
}

// CONSTRUCTOR AND MEMBER FUNCTIONS
MyModel::MyModel()
//:x(nj)
//,y(ni)
//
:blobs(3, 100, false, MyConditionalPrior(), DNest4::PriorType::log_uniform)
{
//    for(size_t i=0; i<ni; ++i)
//        y[i] = y_max - (i + 0.5)*dy;

//    for(size_t j=0; j<nj; ++j)
//        x[j] = x_min + (j + 0.5)*dx;
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

    const auto& _blobs = blobs.get_components();
    double lost_flux = 0.0;
    double R = 1.0; double Rsq = R*R;
    double r, rsq, d, dsq, overlap_area;
    for(size_t b=0; b<_blobs.size(); ++b)
    {
        overlap_area = 0.0;
        r = _blobs[b][3];
        rsq = r*r;
        dsq = _blobs[b][0]*_blobs[b][0] + _blobs[b][1]*_blobs[b][1];
        d = sqrt(dsq);
        if(d < r + R)
        {
            // See Wolfram article at
            // http://mathworld.wolfram.com/Circle-CircleIntersection.html
            overlap_area = rsq*acos((dsq + rsq - Rsq)/(2*d*r))
                                    + Rsq*acos((dsq + Rsq - rsq)/(2*d*R))
                            -0.5*sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R));

            lost_flux += _blobs[b][2]*overlap_area;
        }
        std::cout<<r<<' '<<rsq<<' '<<dsq<<' '<<d<<' '<<overlap_area<<' '<<lost_flux<<std::endl;
    }
    exit(0);

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

} // namespace Obscurity

