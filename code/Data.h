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

    public:
        Data();

        void load(const char* filename);
};

} // namespace Obscurity

#endif

