#pragma once

#include <iostream>
#include "DiscreteFractureNetwork.hpp"

using namespace std;

namespace FractureNetwork {

bool ImportFracture(const string& filePathInput, const string filePathOutput, DiscreteFractureNetwork& fracture, Traces& trace);

bool ReadFracture(const string& filePath, const string& fileName, DiscreteFractureNetwork& fracture);

bool FractureIntersection(const DiscreteFractureNetwork fracture, Traces& trace);

bool FindTraces(const Vector3d s, const Vector3d Point, const DiscreteFractureNetwork& fracture, unsigned int Id1, unsigned int Id2, Traces& trace);

void SaveTraces(double n, double m, Vector3d point, Vector3d s, Traces& trace, unsigned int Id1, unsigned int Id2);

bool PrintOnFile(const string fileName, const string filePath, Traces trace);

}
