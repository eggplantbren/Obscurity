#ifndef DNest4_Template_MyModel
#define DNest4_Template_MyModel

#include "DNest4/code/DNest4.h"
#include <ostream>
#include "Blobby.h"

class MyModel
{
	private:
        // Image size
        static constexpr double x_min = -1.0;
        static constexpr double x_max =  1.0;
        static constexpr double y_min = -1.0;
        static constexpr double y_max =  1.0;
        static constexpr size_t ni = 101;
        static constexpr size_t nj = 101;
        static const double dx, dy;

        // Obscuring blobs
        Blobby blobs;

        // Image of a unit circle obscured by the blobs
        std::vector<std::vector<double>> image;

        // Coordinates of centers of pixels
        std::vector<double> x, y;

        // Calculate the image
        void calculate_image();

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
};

#endif

