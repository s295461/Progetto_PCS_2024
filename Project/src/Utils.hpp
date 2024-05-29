#pragma once

#include <iostream>
#include "DiscreteFractureNetwork.hpp"

using namespace std;

namespace FractureNetwork {

bool ImportFracture(const string& filePath, DiscreteFractureNetwork& fracture, Traces trace);

bool readFracture(const string& filePath, const string& fileName, DiscreteFractureNetwork& fracture);

bool FractureIntersection(const DiscreteFractureNetwork fracture, Traces trace);

bool compare(Vector3d a, Vector3d b);

bool isBetween(Vector3d a, Vector3d b, Vector3d c);

bool findTraces(const Vector3d s, const Vector3d Point, const DiscreteFractureNetwork& fracture, unsigned int Id1, unsigned int Id2, Traces trace);


}
