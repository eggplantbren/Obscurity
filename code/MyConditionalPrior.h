#ifndef Obscurity_MyConditionalPrior
#define Obscurity_MyConditionalPrior

#include "DNest4/code/DNest4.h"

namespace Obscurity
{

class MyConditionalPrior:public DNest4::ConditionalPrior
{
	private:
        // Width of conditional prior for positions
        double sigma;

        // Expected value of conditional prior for "masses" and widths of blobs
        double mu_mass;
        double max_width;

		double perturb_hyperparameters(DNest4::RNG& rng);

	public:
		MyConditionalPrior();

		void from_prior(DNest4::RNG& rng);

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;
		static const int weight_parameter = 1;
};

} // namespace Obscurity

#endif

