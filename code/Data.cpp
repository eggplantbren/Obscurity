#include "Data.h"
#include <fstream>
#include <iostream>
#include <algorithm>

namespace Obscurity
{

Data::Data()
{

}

void Data::load(const char* filename)
{
    std::fstream fin(filename, std::ios::in);
    if(!fin)
    {
        std::cerr<<"# Couldn't open file "<<filename<<'.'<<std::endl;
        return;
    }

    t.clear();
    y.clear();
    sig.clear();

    double temp1, temp2, temp3;
    while(fin>>temp1 && fin>>temp2 && fin>>temp3)
    {
        t.push_back(temp1);
        y.push_back(temp2);
        sig.push_back(temp3);
    }

    t_min = *std::min_element(t.begin(), t.end());
    t_max = *std::max_element(t.begin(), t.end());
    t_range = t_max - t_min;

    std::cout<<"# Loaded "<<t.size()<<" data points from file ";
    std::cout<<filename<<'.'<<std::endl;

    fin.close();
}







} // namespace Obscurity

