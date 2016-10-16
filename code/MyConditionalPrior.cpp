#include "MyConditionalPrior.h"
#include "DNest4/code/DNest4.h"
#include <cmath>
#include <boost/math/distributions/normal.hpp>

namespace Obscurity
{

MyConditionalPrior::MyConditionalPrior()
{

}

void MyConditionalPrior::from_prior(DNest4::RNG& rng)
{
    sigma = rng.rand();

    mu_mass = exp(log(1E-3) + log(1E6)*rng.rand());
    mu_width = exp(log(1E-3) + log(1E6)*rng.rand());
}

double MyConditionalPrior::perturb_hyperparameters(DNest4::RNG& rng)
{
	double logH = 0.0;

    int which = rng.rand_int(3);

    if(which == 0)
    {
        sigma += rng.randh();
        DNest4::wrap(sigma, 0.0, 1.0);
    }
    else if(which == 1)
    {
        mu_mass = log(mu_mass);
        mu_mass += log(1E6)*rng.randh();
        DNest4::wrap(mu_mass, log(1E-3), log(1E3));
        mu_mass = exp(mu_mass);
    }
    else
    {
        mu_width = log(mu_width);
        mu_width += log(1E6)*rng.randh();
        DNest4::wrap(mu_width, log(1E-3), log(1E3));
        mu_width = exp(mu_width);
    }

	return logH;
}

// vec = {xc, yc, mass, width}

double MyConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
    double logp = 0.0;

    if(vec[2] < 0 || vec[3] < 0.999*mu_width || vec[3] > 1.001*mu_width)
        return -1E300;

    logp += -log(2*M_PI*sigma*sigma)
                -0.5*(vec[0]*vec[0] + vec[1]*vec[1])/(sigma*sigma);
    logp += -log(mu_mass) - vec[2]/mu_mass;

	return logp;
}

#include <iostream>
void MyConditionalPrior::from_uniform(std::vector<double>& vec) const
{
    const boost::math::normal standard_normal(0.0, 1.0);
    vec[0] = sigma*quantile(standard_normal, vec[0]);
    vec[1] = sigma*quantile(standard_normal, vec[1]);
    vec[2] = -mu_mass*log(1.0 - vec[2]);
    vec[3] = mu_width*(0.999 + 0.002*vec[3]);
}

void MyConditionalPrior::to_uniform(std::vector<double>& vec) const
{
    const boost::math::normal standard_normal(0.0, 1.0);
    vec[0] = cdf(standard_normal, vec[0]/sigma);
    vec[1] = cdf(standard_normal, vec[1]/sigma);
    vec[2] = 1.0 - exp(-vec[2]/mu_mass);
    vec[3] = (vec[3]/mu_width - 0.999)/0.002;
}

void MyConditionalPrior::print(std::ostream& out) const
{
	out<<sigma<<' '<<mu_mass<<' '<<mu_width<<' ';
}

} // namespace Obscurity

