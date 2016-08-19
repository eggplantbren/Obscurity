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
    static const DNest4::Cauchy c;
    sigma = std::abs(c.generate(rng));

    mu_mass = exp(log(1E-3) + log(1E3)*rng.rand());
    max_width = exp(log(1E-3) + log(1E3)*rng.rand());
}

double MyConditionalPrior::perturb_hyperparameters(DNest4::RNG& rng)
{
	double logH = 0.0;

    int which = rng.rand_int(3);
    static const DNest4::Cauchy c;

    if(which == 0)
    {
        logH += c.perturb(sigma, rng);
        sigma = std::abs(sigma);
    }
    else if(which == 1)
    {
        mu_mass = log(mu_mass);
        mu_mass += log(1E3)*rng.randh();
        DNest4::wrap(mu_mass, log(1E-3), log(1.0));
        mu_mass = exp(mu_mass);
    }
    else
    {
        max_width = log(max_width);
        max_width += log(1E3)*rng.randh();
        DNest4::wrap(max_width, log(1E-3), log(1.0));
        max_width = exp(max_width);
    }

	return logH;
}

// vec = {xc, yc, mass, width}

double MyConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
    double logp = 0.0;

    if(vec[2] < 0 || vec[3] < 0 || vec[3] > max_width)
        return -1E300;

    logp += -log(2*M_PI*sigma*sigma)
                -0.5*(vec[0]*vec[0] + vec[1]*vec[1])/(sigma*sigma);
    logp += -log(mu_mass) - vec[2]/mu_mass;
    logp += -log(max_width);    

	return logp;
}

void MyConditionalPrior::from_uniform(std::vector<double>& vec) const
{
    static const boost::math::normal standard_normal(0.0, 1.0);
    vec[0] = sigma*quantile(standard_normal, vec[0]);
    vec[1] = sigma*quantile(standard_normal, vec[1]);
    vec[2] = -mu_mass*log(1.0 - vec[2]);
    vec[3] = max_width*vec[3];
}

void MyConditionalPrior::to_uniform(std::vector<double>& vec) const
{
    static const boost::math::normal standard_normal(0.0, 1.0);
    vec[0] = cdf(standard_normal, vec[0]/sigma);
    vec[1] = cdf(standard_normal, vec[1]/sigma);
    vec[2] = 1.0 - exp(-vec[2]/mu_mass);
    vec[3] = vec[3]/max_width;
}

void MyConditionalPrior::print(std::ostream& out) const
{
	out<<' ';
}

} // namespace Obscurity

