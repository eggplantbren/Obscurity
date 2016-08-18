#ifndef Obscurity_Data
#define Obscurity_Data

#include <vector>

namespace Obscurity
{

class Data
{
    private:
        std::vector<double> t;
        std::vector<double> y;
        std::vector<double> sig;

        double t_min, t_max, t_range;

    public:
        Data();

        void load(const char* filename);

        // Getters
        const std::vector<double>& get_t() const
        { return t; }
        const std::vector<double>& get_y() const
        { return y; }
        const std::vector<double>& get_sig() const
        { return sig; }

        double get_t_min() const
        { return t_min; }
        double get_t_max() const
        { return t_max; }
        double get_t_range() const
        { return t_range; }
};

} // namespace Obscurity

#endif

