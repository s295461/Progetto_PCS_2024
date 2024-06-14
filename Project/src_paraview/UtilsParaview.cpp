#include <UtilsParaview.hpp>

#include <fstream>

namespace GeometryLibrary{
//*********************************************************
vector<vector<vector<unsigned int>>> Polygons::TriangulatePolygons()
{
    const unsigned int numPolygons = listVertices.size();
    vector<vector<vector<unsigned int>>> triangleList(numPolygons);

    for(unsigned int p = 0; p < numPolygons; p++)
    {
        const unsigned int numPolygonVertices = listVertices[p].size();

        for (unsigned int v = 0; v < numPolygonVertices; v++)
        {
            const unsigned int nextVertex = listVertices[p][(v + 1) % numPolygonVertices];
            const unsigned int nextNextVertex = listVertices[p][(v + 2) % numPolygonVertices];

            if ((v + 2) % numPolygonVertices == 0)
                break;

            vector<unsigned int> triangle_vertices = {listVertices[p][0], nextVertex, nextNextVertex};

            triangleList[p].push_back(triangle_vertices);
        }
    }
    return triangleList;
}
//*********************************************************
vector<double> Polygons::computePolygonsArea()
{


    vector<vector<vector<unsigned int>>> triangleList = TriangulatePolygons();

    const unsigned int numPolygons = listVertices.size();
    vector<double> area(numPolygons, 0.0);
    for(unsigned int p = 0; p < numPolygons; p++)
    {
        for(unsigned int t = 0; t < triangleList.size(); t++)
        {
            Matrix3d points;
            points << VerticesCoordinates.col(triangleList[p][t][0]),  VerticesCoordinates.col(triangleList[p][t][1]), VerticesCoordinates.col(triangleList[p][t][2]);

            Triangle triangle(points);

            area[p] += triangle.computeArea();
        }
    }

    return area;
}
//*********************************************************
void Polygons::GedimInterface(vector<vector<unsigned int>>& triangles,
                              VectorXi& materials)
{
    const unsigned int numPolygons = listVertices.size();
    vector<vector<vector<unsigned int>>> triangleList = TriangulatePolygons();

    unsigned int numTotalTriangles = 0;
    for(unsigned int p = 0; p < numPolygons; p++)
        numTotalTriangles += triangleList[p].size();

    triangles.reserve(numTotalTriangles);
    materials = VectorXi::Zero(numTotalTriangles);

    unsigned int count = 0;
    for(unsigned int p = 0; p < numPolygons; p++)
    {
        for(unsigned int t = 0; t < triangleList[p].size(); t++)
        {
            triangles.push_back(triangleList[p][t]);
            materials(count) = p;
            count++;
        }
    }
}
//*********************************************************
void importPolygonsList(const string& filepath,
                        Polygons &polygons)
{
    ifstream file(filepath);

    if (file.fail())
        throw runtime_error("cannot open file");

    string line;

    // Read number of fractures
    getline(file, line); // header
    getline(file, line);
    unsigned int numpolygons;
    istringstream convertNP(line);
    convertNP >> numpolygons;

    polygons.listVertices.resize(numpolygons);

    // Loop through each fracture
    for (unsigned int c = 0; c < numpolygons; c++) {
        getline(file, line); // header "FractureId; NumVertices"
        getline(file, line);

        unsigned int id, numvertices;
        char separator;
        istringstream convertNV(line);
        convertNV >> id >> separator >> numvertices;

        polygons.listVertices[c].resize(numvertices);

        getline(file, line); // header "Vertices"

        MatrixXd vertices(3, numvertices);

        for (unsigned int i = 0; i < 3; i++) {
            getline(file, line);
            istringstream convertV(line);
            for (unsigned int v = 0; v < numvertices; v++) {
                double coord;
                convertV >> coord >> separator;
                vertices(i, v) = coord;
            }
        }

        polygons.VerticesCoordinates.conservativeResize(3, polygons.VerticesCoordinates.cols() + numvertices);
        polygons.VerticesCoordinates.rightCols(numvertices) = vertices;

        // Update listVertices with the indices of the vertices for this fracture
        for (unsigned int v = 0; v < numvertices; v++) {
            polygons.listVertices[c][v] = polygons.VerticesCoordinates.cols() - numvertices + v;
        }
    }

    file.close();
}
//********************************************************
}
