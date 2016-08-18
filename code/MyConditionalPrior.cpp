#include "MyConditionalPrior.h"
#include "DNest4/code/DNest4.h"
#include <cmath>

namespace Obscurity
{

MyConditionalPrior::MyConditionalPrior()
{

}

void MyConditionalPrior::from_prior(DNest4::RNG& rng)
{
    DNest4::Cauchy c;
    sigma = std::abs(c.generate(rng));

    mu_mass = exp(log(1E-3) + log(1E3)*rng.rand());
    mu_width = exp(log(1E-3) + log(1E3)*rng.rand());
}

double MyConditionalPrior::perturb_hyperparameters(DNest4::RNG& rng)
{
	double logH = 0.0;



	return logH;
}

// vec = {xc, yc, mass, width}

double MyConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
	return 0.;
}

void MyConditionalPrior::from_uniform(std::vector<double>& vec) const
{
    vec[0] = 
}

void MyConditionalPrior::to_uniform(std::vector<double>& vec) const
{

}

void MyConditionalPrior::print(std::ostream& out) const
{
	out<<' ';
}

} // namespace Obscurity

