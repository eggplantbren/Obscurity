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

	double rsq, widthsq, tausq;
    double C = 2.0/M_PI;
	for(size_t i=0; i<components.size(); i++)
	{
		rsq = pow(x - components[i][0], 2)
				+ pow(y - components[i][1], 2);

		widthsq = components[i][3]*components[i][3];

		if(rsq < widthsq)
        {
            tausq = 1.0/widthsq;
			f += C*components[i][2]*(1.0 - rsq*tausq)*tausq;
        }
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

