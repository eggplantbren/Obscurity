#ifndef Obscurity_MyModel_Pixels
#define Obscurity_MyModel_Pixels

#include "DNest4/code/DNest4.h"
#include <ostream>
#include <armadillo>
#include "Data.h"
#include "MyConditionalPrior.h"

namespace Obscurity
{

class MyModel_Pixels
{
	private:
        // Coordinates for pre-whitened parameterisation of obscuring map
        arma::mat n;

        // Scale sd for creating obscurer map
        double sd;

        // Obscuring map
        arma::mat obscurer_map;
        arma::mat convolved;

        // Initial positional offset of blobs, and their crossing timescale
        double x0, timescale;

        // Calculate the obscurer map
        void calculate_obscurer_map();

        // Calculate the total flux
        double calculate_total_flux(double time) const;

	public:
		// Constructor only gives size of params
		MyModel_Pixels();

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


        /* STATIC STUFF BEGINS HERE */
    private:
        static Data data;

        // Image size parameters
        static constexpr size_t ni = 101;
        static constexpr size_t nj = 201;
        static constexpr double x_min = -4.0;
        static constexpr double x_max =  4.0;
        static constexpr double y_min = -2.0;
        static constexpr double y_max =  2.0;
        static constexpr double dx = (x_max - x_min)/nj;
        static constexpr double dy = (y_max - y_min)/ni;
        static constexpr double LL = sqrt(dx*dy);

        // Star image
        static arma::mat star;
        static arma::cx_mat fft_of_star;

        // Coordinates of pixel centers
        static std::vector<double> x, y;

    public:
        static void initialise();
        static void load_data(const char* filename);
};

} // namespace Obscurity

#endif

