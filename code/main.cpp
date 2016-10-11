#include <iostream>
#include "DNest4/code/DNest4.h"
#include "MyModel.h"
#include "MyModel_Pixels.h"

// Which model am I using?
typedef Obscurity::MyModel_Pixels TheModel;

int main(int argc, char** argv)
{
    TheModel::load_data("data.txt");
    TheModel::initialise();

    DNest4::start<TheModel>(argc, argv);
    return 0;
}

