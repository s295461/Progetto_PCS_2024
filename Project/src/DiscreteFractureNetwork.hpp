#pragma once

#include <iostream>
#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace FractureNetwork {

struct DiscreteFractureNetwork
{
    unsigned int numFracture = 0;
    vector<unsigned int> fractureID = {};
    vector<unsigned int> NumVertices = {};
    vector<MatrixXd> vertices;
};

// struct Vec3d
// {
//     double x, y, z;
// };

struct BoundingBox {
    Vector3d min;
    Vector3d max;
};
struct FractureBBox{
    BoundingBox bbox;
    unsigned int fractureID;
};
}
