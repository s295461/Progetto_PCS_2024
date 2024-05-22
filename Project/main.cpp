#include <iostream>
#include <vector>
#include "Eigen/Eigen"
#include "Utils.hpp"
#include "DiscreteFractureNetwork.hpp"


using namespace std;
using namespace Eigen;
using namespace FractureNetwork;

// struct Point {
//     double x, y, z;
// };

// // Funzione che calcola il prodotto vettoriale
// Point crossProduct(const Point& a, const Point& b)
// {
//     Point result;
//     result.x = a.y * b.z - a.z * b.y;
//     result.y = a.z * b.x - a.x * b.z;
//     result.z = a.x * b.y - a.y * b.x;
//     return result;

// }

// // Funzione che calcola il prodotto scalare di due vettori
// double dotProduct(const Point& a, const Point& b)
// {
//     return a.x * b.x + a.y * b.y + a.z * b.z;
// }

// // Funzione che calcola la norma di un vettore
// double norm(const Point& a)
// {
//     return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
// }

// // Funzione per calcolare il punto di intersezione tra due segmenti
// Point findIntersection(const Point& p1, const Point& p2, const Point& q1, const Point& q2)
// {
//     Point direction1 = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
//     Point direction2 = {q2.x - q1.x, q2.y - q1.y, q2.z - q1.z};

//     Point normal;
//     normal.x = crossProduct(direction1, direction2).x / (norm(direction1) * norm(direction2));
//     normal.y = crossProduct(direction1, direction2).y / (norm(direction1) * norm(direction2));
//     normal.z = crossProduct(direction1, direction2).z / (norm(direction1) * norm(direction2));

//     if(abs(normal.x) == 0 && abs(normal.y) == 0 && abs(normal.z) == 0)
//         return {INFINITY, INFINITY, INFINITY};

//     Point p1q1 = {q1.x - p1.x, q1.y - p1.y, q1.z - p1.z};

//     double t = dotProduct(p1q1, normal) / dotProduct(direction1, normal);

//     Point intersection;
//     intersection.x = p1.x + direction1.x * t;
//     intersection.y = p1.y + direction1.y * t;
//     intersection.z = p1.z + direction1.z * t;

//     return intersection;
// }

// //Funzione che trova i punti di intersezione tra due poligoni
// vector<Point> findPolygonIntersections(const vector<Point>& poly1, const vector<Point> poly2)
// {
//     vector<Point> intersections;

//     //Iteriamo su tutti i lati di entrambi i poligoni
//     for(unsigned int i = 0; i < poly1.size(); ++i)
//     {
//         unsigned int next_i = (i + 1) % poly1.size();
//         for(unsigned int j = 0; j < poly2.size(); ++j)
//         {
//             unsigned int next_j = (j + 1) % poly2.size();
//             // Verifico l'intersezione tra due segmenti
//             Point intersection = findIntersection(poly1[i], poly1[next_i], poly2[j], poly2[next_j]);

//             // Verifico se il punto di intersezione si trova all'interno di entrambi i segmenti
//             if(intersection.x != INFINITY && intersection.y != INFINITY && intersection.z != INFINITY&&
//                 ((intersection.x >= min(poly1[i].x, poly1[next_i].x)&&
//                 intersection.x <= max(poly1[i].x, poly1[next_i].x) &&
//                 intersection.y >= min(poly1[i].y, poly1[next_i].y)&&
//                 intersection.y <= max(poly1[i].y, poly1[next_i].y) &&
//                 intersection.z >= min(poly1[i].z, poly1[next_i].z) &&
//                 intersection.z <= max(poly1[i].z, poly1[next_i].z)) ||
//                 (intersection.x >= min(poly2[j].x, poly2[next_j].x)&&
//                 intersection.x <= max(poly2[j].x, poly2[next_j].x) &&
//                 intersection.y >= min(poly2[j].y, poly2[next_j].y)&&
//                 intersection.y <= max(poly2[j].y, poly2[next_j].y) &&
//                 intersection.z >= min(poly2[j].z, poly2[next_j].z) &&
//                 intersection.z <= max(poly2[j].z, poly2[next_j].z))))
//             {
//                 intersections.push_back(intersection);
//             }

//         }
//     }
//     return intersections;
// }


int main()
{
    DiscreteFractureNetwork fracture;
    string filePath = "DFN";
    if(!ImportFracture(filePath, fracture))
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

// branch sara
