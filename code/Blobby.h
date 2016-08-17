#ifndef Lensing2_Blobby
#define Lensing2_Blobby

#include "DNest4/code/DNest4.h"

class Blobby
{
	private:
		DNest4::RJObject<DNest4::BasicCircular> blobs;

	public:
		// Pass in image scale and flux scale
		Blobby(double x_min, double x_max,
					double y_min, double y_max);

		// Required methods
		double evaluate(double x, double y, bool update=false) const;
		void from_prior(DNest4::RNG& rng);
		double perturb(DNest4::RNG& rng);
		void print(std::ostream& out) const;

        // Getter
		const DNest4::RJObject<DNest4::BasicCircular>& get_blobs() const
        { return blobs; }
};

#endif

