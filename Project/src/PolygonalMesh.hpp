#pragma once

#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


struct PolygonalMesh {

    // Cell0D
    unsigned int numCell0D = 0;
    vector<unsigned int> cellId0D = {};
    vector<Vector3d> coordinates0D = {};


    // Cell1D
    unsigned int numCell1D = 0;
    vector<unsigned int> cellId1D = {};
    vector<Vector2i> verticesId1D = {};


    // Cell2D
    unsigned int numCell2D = 0;
    vector<unsigned int> cellId2D = {};
    vector<unsigned int> numVertices2D = {};
    vector<unsigned int> numEdges2D = {};
    vector<vector<unsigned int>> verticesId2D = {};
    vector<vector<unsigned int>> edgesId2D = {};
};



