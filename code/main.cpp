#include <iostream>
#include "DNest4/code/DNest4.h"
#include "MyModel.h"

// Which model am I using?
typedef Obscurity::MyModel TheModel;

int main(int argc, char** argv)
{
    TheModel::load_data("easy_data.txt");
    TheModel::initialise();

    DNest4::start<TheModel>(argc, argv);
    return 0;
}

