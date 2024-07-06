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
using namespace GeometryLibrary;

int main(){
    string filepath = "C:/Users/davep/OneDrive/Desktop/Progetto_PCS_2024/Project/DFN/FR3_data.txt";
    string filepath2 = "C:/Users/davep/OneDrive/Desktop/Progetto_PCS_2024/Project/Result/FR82_traces.txt";

    GeometryLibrary::Polygons polygons;
    GeometryLibrary::importPolygonsList(filepath, polygons);
    Eigen::VectorXi materials;

    Gedim::UCDUtilities exporter;
    std::vector<std::vector<unsigned int>> triangles;

    // Importa i poligoni
    polygons.GedimInterface(triangles, materials);
    exporter.ExportPolygons("C:/Users/davep/OneDrive/Desktop/Progetto_PCS_2024/Project/Export_Paraview/FR3_paraview.inp",
                            polygons.VerticesCoordinates,
                            triangles,
                            {},
                            {},
                            materials);

    MatrixXd points(3,0);
    MatrixXi index_edges(2,0);
    VectorXi materials2;

    // Importa i segmenti
    importSegments(filepath2, points, index_edges, materials2);

    exporter.ExportSegments("C:/Users/davep/OneDrive/Desktop/Progetto_PCS_2024/Project/Export_Paraview/FR82_traces_paraview.inp",
                            points,
                            index_edges,
                            {},
                            {},
                            materials2);


    return 0;
}
