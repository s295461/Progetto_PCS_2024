#pragma once

#include <iostream>
#include "DiscreteFractureNetwork.hpp"

using namespace std;

namespace FractureNetwork {

bool ImportFracture(const string& filePath, DiscreteFractureNetwork& fracture);

bool readFracture(const string& filePath, const string& fileName, DiscreteFractureNetwork& fracture);

// Vec3d crossProduct(const Vec3d& a, const Vec3d& b);

// double norm2(const Vec3d& a);

// bool findTraces(const Vec3d s, const Vec3d Point, const DiscreteFractureNetwork& fracture, unsigned int Id1, unsigned int Id2);

bool PlaneIntersection(const DiscreteFractureNetwork fracture);



}
