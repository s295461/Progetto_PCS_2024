#pragma once

#include <iostream>
#include "DiscreteFractureNetwork.hpp"

using namespace std;

namespace FractureNetwork {

bool ImportFracture(const string& filePath, DiscreteFractureNetwork fracture);

bool readFracture(const string& filePath, const string& fileName, DiscreteFractureNetwork fracture);


}
