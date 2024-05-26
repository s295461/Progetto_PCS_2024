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
    DiscreteFractureNetwork fracture;
    Traces trace;
    string filePath = "DFN";
    if(!ImportFracture(filePath, fracture, trace))
        return 1;



    // PlaneIntersection(fracture);




    /// PROVE
    // Vector3i P = {0, 0, 0};
    // Vector3i Q = {1, 0, 0};
    // Vector3i R = {1, 1, 0};
    // Vector3i S = {0, 1, 0};

    // int a = 1;
    // int a0 = 1;
    // int a1 = 1;

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


    // // Poligoni:
    // vector<Point> poly1 = {{0, 0, 0}, {4, 0, 0}, {4, 4, 0}, {0, 4, 0}};
    // vector<Point> poly2 = {{2, 1, 1}, {4, 1, 1}, {4, 1, -2}, {2, 1, -2}};

    // // Trova i punti di intersezione:
    // vector<Point> intersections = findPolygonIntersections(poly1, poly2);

    // // Stampo i punti di intersezione
    // cout << "Intersection points:" << endl;
    // for(const auto& point : intersections)
    //     cout << "(" << point.x << ", " << point.y << ", " << point.z << ")" << endl;

  return 0;
}
