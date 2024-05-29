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

struct Traces
{
    unsigned int numTraces = 0;
    vector<unsigned int> traceId = {};
    vector<Vector2i> fractureId = {};
    vector<MatrixXd> coordinates = {};
};


}


//metterci dentro ID --> ottengo la matrice dei vertici
