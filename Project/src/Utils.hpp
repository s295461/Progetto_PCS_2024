#pragma once

#include "DiscreteFractureNetwork.hpp"
#include "PolygonalMesh.hpp"

using namespace std;

namespace FractureNetwork {

bool ImportFracture(const string fileNameInput, const string fileNameOutput, const string fileNameOutputReordered, const string filePathInput, const string filePathOutput, DiscreteFractureNetwork& fracture, Traces& trace);

void clearDiscreteFractureNetwork(DiscreteFractureNetwork& fracture);

void clearTraces(Traces& trace);


bool ReadFracture(const string& filePath, const string& fileName, DiscreteFractureNetwork& fracture);


BoundingBox BBox3D(const MatrixXd& vertices);

bool FractureIntersection(const DiscreteFractureNetwork fracture, Traces& trace);


void SaveTraces(double n, double m, Vector3d point, Vector3d s, Traces& trace, unsigned int Id1, unsigned int Id2);

bool PrintOnFile(const string fileName, const string filePath, Traces trace);

bool TraceReorder(DiscreteFractureNetwork& fracture, Traces& trace);

double PointDistance(Vector3d P, Vector3d Q);

bool reordering(vector<unsigned int>& idTraces, vector<double>& length);

bool printTraces(const string fileName, const string filePath, Traces trace, DiscreteFractureNetwork fracture);

bool FindTraces(const Vector3d s, const Vector3d Point, const DiscreteFractureNetwork& fracture, unsigned int Id1, unsigned int Id2, Traces& trace);

void SaveTraces(double n, double m, Vector3d point, Vector3d s, Traces& trace, unsigned int Id1, unsigned int Id2);

bool PrintOnFile(const string fileName, const string filePath, Traces trace);


bool TraceReorder(DiscreteFractureNetwork& fracture, Traces& trace);

double PointDistance(Vector3d P, Vector3d Q);

bool reordering(vector<unsigned int>& idTraces, vector<double>& length);

bool printTraces(const string fileName, const string filePath, Traces trace, DiscreteFractureNetwork fracture);

void rewriteData(const string& inputFileName, const string& outputFileName);


}

namespace PolygonalMesh {

bool fractureCut(FractureNetwork::DiscreteFractureNetwork& fracture, FractureNetwork::Traces& trace, Cell0D& Cell0D, Cell1D& Cell1D, Cell2D& Cell2D);

bool createSubfracture(vector<Vector3d> subfracture, vector<Vector3d> cuttingTrace, vector<vector<Vector3d>>& subfractureVertices);

vector<Vector3d> extendTraces(MatrixXd fractureVertices, MatrixXd verticesTrace);

// vector<Vector3d> intersectTraces(FractureNetwork::Traces trace, unsigned int idTrace, unsigned int idFracture);

vector<tuple<Vector3d, unsigned int, unsigned int>> intersectTraces(FractureNetwork::Traces trace, unsigned int idFracture);

}
