#include <iostream>
#include <vector>
#include "Eigen/Eigen"
#include "Utils.hpp"
#include "DiscreteFractureNetwork.hpp"


using namespace std;
using namespace Eigen;
using namespace FractureNetwork;

int main()
{
    DiscreteFractureNetwork fracture;
    Traces trace;
    string filePathInput = "DFN";
    string filePathOutput = "Result";

    if(!ImportFracture(filePathInput, filePathOutput, fracture, trace))
        return 1;

  return 0;
}
