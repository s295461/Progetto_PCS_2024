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
void importPolygonsList(const string& filepath, Polygons &polygons)
{
    // Apre il file
    ifstream file(filepath);

    // Controlla se l'apertura del file è fallita
    if (file.fail())
        throw runtime_error("cannot open file");

    string line;

    // Legge e ignora l'header del file
    getline(file, line);
    getline(file, line);

    // Legge il numero di poligoni
    unsigned int numpolygons;
    istringstream convertNP(line);
    convertNP >> numpolygons;

    // Ridimensiona la lista dei vertici per contenere il numero di poligoni
    polygons.listVertices.resize(numpolygons);

    // Itera su ciascun poligono
    for (unsigned int c = 0; c < numpolygons; c++) {
        // Legge e ignora l'header "FractureId; NumVertices"
        getline(file, line);
        getline(file, line);

        unsigned int id, numvertices;
        char separator;
        istringstream convertNV(line);
        convertNV >> id >> separator >> numvertices;

        // Ridimensiona la lista dei vertici per il poligono corrente
        polygons.listVertices[c].resize(numvertices);

        // Legge e ignora l'header "Vertices"
        getline(file, line);

        // Matrice per contenere le coordinate dei vertici
        MatrixXd vertices(3, numvertices);

        // Legge le coordinate dei vertici
        for (unsigned int i = 0; i < 3; i++) {
            getline(file, line);
            istringstream convertV(line);
            for (unsigned int v = 0; v < numvertices; v++) {
                double coord;
                convertV >> coord >> separator;
                vertices(i, v) = coord;
            }
        }

        // Aggiunge le nuove coordinate alla matrice esistente dei vertici
        polygons.VerticesCoordinates.conservativeResize(3, polygons.VerticesCoordinates.cols() + numvertices);
        polygons.VerticesCoordinates.rightCols(numvertices) = vertices;

        // Aggiorna la lista dei vertici per il poligono corrente
        for (unsigned int v = 0; v < numvertices; v++) {
            polygons.listVertices[c][v] = polygons.VerticesCoordinates.cols() - numvertices + v;
        }
    }

    file.close();
}
//********************************************************
void importSegments(const string& filePath, MatrixXd& points, MatrixXi& index_edges, VectorXi& materials) {
    // Apre il file
    ifstream file(filePath);

    // Controlla se il file è stato aperto correttamente
    if (!file.is_open()) {
        cerr << "Impossibile aprire il file " << filePath << endl;
        return;
    }

    string line;
    int lineNumber = 0;

    // Ignora le prime due righe del file
    getline(file, line);
    lineNumber++;
    getline(file, line);
    lineNumber++;

    // Legge i segmenti dal file
    while (getline(file, line)) {
        lineNumber++;
        // Salta le righe vuote o i commenti
        if (line.empty() || line[0] == '#')
            continue;

        istringstream iss(line);
        int traceId, fractureId1, fractureId2;
        double x1, y1, z1, x2, y2, z2;
        char delimiter;

        // Legge i dati del segmento
        iss >> traceId >> delimiter >> fractureId1 >> delimiter >> fractureId2
            >> delimiter >> x1 >> delimiter >> y1 >> delimiter >> z1
            >> delimiter >> x2 >> delimiter >> y2 >> delimiter >> z2;

        // Aggiunge i punti alla matrice dei punti
        points.conservativeResize(3, points.cols() + 2);
        points.col(points.cols() - 2) << x1, y1, z1;
        points.col(points.cols() - 1) << x2, y2, z2;

        // Aggiunge gli indici degli edge alla matrice degli indici
        index_edges.conservativeResize(2, index_edges.cols() + 1);
        index_edges.col(index_edges.cols() - 1) << points.cols() - 2, points.cols() - 1;

        // Aggiunge il materiale al vettore dei materiali
        materials.conservativeResize(materials.size() + 1);
        materials(materials.size() - 1) = traceId;
    }

    file.close();
}
}

