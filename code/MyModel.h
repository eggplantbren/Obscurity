#ifndef Obscurity_MyModel
#define Obscurity_MyModel

#include "DNest4/code/DNest4.h"
#include <ostream>
#include <tuple>
#include "Data.h"
#include "Blobby.h"

namespace Obscurity
{

class MyModel
{
	private:
        static Data data;

        // Image size
        static constexpr double x_min = -1.0;
        static constexpr double x_max =  1.0;
        static constexpr double y_min = -1.0;
        static constexpr double y_max =  1.0;
        static constexpr size_t ni = 201;
        static constexpr size_t nj = 201;
        static const double dx, dy;

        // Obscuring blobs
        Blobby blobs;

        // Initial positional offset of blobs, and their crossing timescale
        double x0, timescale;

        // Calculate the total flux
        double calculate_total_flux(double time) const;

	public:
		// Constructor only gives size of params
		MyModel();

		// Generate the point from the prior
		void from_prior(DNest4::RNG& rng);

		// Metropolis-Hastings proposals
		double perturb(DNest4::RNG& rng);

		// Likelihood function
		double log_likelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;

        static void load_data(const char* filename);
};

} // namespace Obscurity

#endif

