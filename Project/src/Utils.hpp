#pragma once

#include "DiscreteFractureNetwork.hpp"
#include "PolygonalMesh.hpp"

using namespace std;


bool ImportFracture(const string fileNameInput, const string fileNameOutput, const string fileNameOutputReordered,
                    const string filePathInput, const string filePathOutput, DiscreteFractureNetwork& fracture, Traces& trace);

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

bool reordering(vector<unsigned int>& idTraces, vector<double>& length);

bool printTraces(const string fileName, const string filePath, Traces trace, DiscreteFractureNetwork fracture);


// *****************************************************************************************************************************************************


bool fractureCut(DiscreteFractureNetwork& fracture, Traces& trace);

bool createSubfracture(vector<Vector3d> subfracture, vector<Vector3d> cuttingTrace, vector<vector<Vector3d>>& subfractureVertices);

vector<Vector3d> extendTraces(vector<Vector3d> subFractureVertices, vector<Vector3d> subTraceVertices);

bool createMesh(vector<pair<vector<Vector3d>, vector<pair<vector<Vector3d>, unsigned int>>>> subFracture, PolygonalMesh& mesh);
