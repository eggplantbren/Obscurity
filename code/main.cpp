#include <iostream>
#include "DNest4/code/DNest4.h"
#include "MyModel.h"

int main(int argc, char** argv)
{
    Obscurity::MyModel::load_data("data.txt");
    DNest4::start<Obscurity::MyModel>(argc, argv);
    return 0;
}

