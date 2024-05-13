#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace FractureNetwork {

struct DiscreteFractureNetwork
{
    unsigned int numFracture = 0;
    vector<unsigned int> fractureID = {};
    vector<Vector2d> coordinates = {};



};




}
