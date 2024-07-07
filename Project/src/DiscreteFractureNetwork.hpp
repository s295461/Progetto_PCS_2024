#pragma once

#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

// Struttura per memorizzare le fratture
struct DiscreteFractureNetwork
{
    unsigned int numFracture = 0;
    vector<unsigned int> fractureID = {};
    vector<unsigned int> NumVertices = {};
    vector<MatrixXd> vertices = {};
};

// Struttura per memorizzare le tracce
struct Traces
{
    unsigned int numTraces = 0;
    vector<unsigned int> traceId = {};
    vector<Vector2i> fractureId = {};
    vector<MatrixXd> coordinates = {};

    vector<double> length = {};
    vector<vector<tuple<unsigned int, bool, double>>> traceReordered;
};



