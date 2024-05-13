#include <iostream>
#include <vector>
#include "Eigen/Eigen"
#include "Utils.hpp"
#include "DiscreteFractureNetwork.hpp"


using namespace std;
using namespace Eigen;
using namespace FractureNetwork;


int main()
{
    string filePath = "DFN";
    ImportFracture(filePath);




    /// PROVE
    // Vector3i P = {5, 2, 3};
    // Vector3i Q = {4, 7, 1};
    // Vector3i R = {2, 6, 8};

    // int a = 2;
    // int a0 = 3;
    // int a1 = 6;

    // Vector3i x = P + a*(Q - P);
    // Vector3i x1 = a*Q + (1-a)*P;

    // for(unsigned int i = 0; i<3; i++)
    //     cout << x[i] << endl;

    // for(unsigned int i = 0; i<3; i++)
    //     cout << x1[i] << endl;


    // Vector3i y = a0*P + a1*Q + (1-a0-a1)*R;
    // for(unsigned int i = 0; i<3; i++)
    //     cout << y[i] << endl;

    // Vector3i u = P - R;
    // Vector3i v = Q - R;
    // Vector3d n;
  return 0;
}
