#include "Blobby.h"
#include "DNest4/code/RNG.h"

namespace Obscurity
{

Blobby::Blobby(double x_min, double x_max,
					double y_min, double y_max)
:blobs(4, 100, false,
	DNest4::BasicCircular(x_min, x_max, y_min, y_max, 1E-3, 1E3),
            DNest4::PriorType::log_uniform)
{

}

double Blobby::evaluate(double x, double y, bool update) const
{
	double f = 0.0;

	const auto& components = (update)?
                             (blobs.get_added()):
                             (blobs.get_components());

	double rsq, widthsq;
	for(size_t i=0; i<components.size(); i++)
	{
		rsq = pow(x - components[i][0], 2)
				+ pow(y - components[i][1], 2);
		widthsq = pow(components[i][3], 2);

		if(rsq < widthsq)
			f += components[i][2]*
				2./M_PI*(1. - rsq/widthsq)/widthsq;
	}

	return f;
}

void Blobby::from_prior(DNest4::RNG& rng)
{
	blobs.from_prior(rng);
}

double Blobby::perturb(DNest4::RNG& rng)
{
	double logH = 0.;

	logH += blobs.perturb(rng);

	return logH;
}

void Blobby::print(std::ostream& out) const
{
	blobs.print(out);
}

} // namespace Obscurity

