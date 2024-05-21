#pragma once

#include <iostream>
#include "DiscreteFractureNetwork.hpp"

using namespace std;

namespace FractureNetwork {

bool ImportFracture(const string& filePath, DiscreteFractureNetwork& fracture);

bool readFracture(const string& filePath, const string& fileName, DiscreteFractureNetwork& fracture);

Vec3d crossProduct(const Vec3d& a, const Vec3d& b);

double norm2(const Vec3d& a);

bool PlaneIntersection(const DiscreteFractureNetwork fracture);

}
