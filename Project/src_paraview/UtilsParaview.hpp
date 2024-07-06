#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "C:\Users\davep\OneDrive\Desktop\Progetto_PCS_2024\Project\src_paraview\UCDUtilities.hpp"
#include <DiscreteFractureNetwork.hpp>
#include <Eigen/Eigen>

using namespace Eigen;
using namespace std;

namespace GeometryLibrary {

struct Triangle{

    Matrix3d Vertices;

    Triangle(const MatrixXd& Vertices): Vertices(Vertices) {}

    double computeArea()
    {
        double area = 0.0;
        for(unsigned int i = 0; i < 3; i++)
            area += Vertices(0, i) * Vertices(1, (i + 1) % 3) - Vertices(1, i) * Vertices(0, (i + 1) % 3);

        area *= 0.5;

        return area;
    }
};

struct Polygons{

    MatrixXd VerticesCoordinates;
    vector<vector<unsigned int>> listVertices;

    Polygons() = default;
    Polygons(const MatrixXd& VerticesCoordinates,
             const vector<vector<unsigned int>>& listVertices):
        VerticesCoordinates(VerticesCoordinates),
        listVertices(listVertices)
    {}

    vector<vector<vector<unsigned int>>> TriangulatePolygons();
    vector<double> computePolygonsArea();

    void GedimInterface(vector<vector<unsigned int>>& triangles,
                        VectorXi& materials);

};

struct Segment {
    int id1;
    int id2;
    int material;
};

void importPolygonsList(const string& filepath,
                        Polygons& polygons);
void importSegments(const std::string& filePath, MatrixXd& points, MatrixXi& index_edges, VectorXi& materials);
}


