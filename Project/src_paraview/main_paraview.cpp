#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <UCDUtilities.hpp>
#include <UtilsParaview.hpp>
#include <DiscreteFractureNetwork.hpp>


using namespace Eigen;
using namespace std;
using namespace Gedim;

int main(){
    string filepath = "C:/Users/davep/OneDrive/Desktop/Progetto_PCS_2024/Project/DFN/FR200_data.txt";
    GeometryLibrary::Polygons polygons;

    GeometryLibrary::importPolygonsList(filepath, polygons);

    Gedim::UCDUtilities exporter;
    std::vector<std::vector<unsigned int>> triangles;
    Eigen::VectorXi materials;
    polygons.GedimInterface(triangles, materials);
    exporter.ExportPolygons("C:/Users/davep/OneDrive/Desktop/Progetto_PCS_2024/Project/Export_Paraview/FR200_paraview.inp",
                            polygons.VerticesCoordinates,
                            triangles,
                            {},
                            {},
                            materials);

    return 0;
}
