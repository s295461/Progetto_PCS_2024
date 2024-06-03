#pragma once

#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

// namespace PolygonalMesh {

struct Cell0D
{
    unsigned int numCell = 0;
    vector<unsigned int> cellId = {};
    vector<Vector3d> coordinates = {};
};

struct Cell1D
{
    unsigned int numCell = 0;
    vector<unsigned int> cellId = {};
    vector<Vector2i> verticesId = {};
};

struct Cell2D
{
    unsigned int numCell = 0;
    vector<unsigned int> numVertices = {};
    vector<unsigned int> numEdges = {};
    vector<list<unsigned int>> verticesId = {};
    vector<list<unsigned int>> edgesId = {};
};

// }
