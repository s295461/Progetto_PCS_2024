#ifndef TEST_H
#define TEST_H

#include <iostream>
#include "Eigen/Eigen"
#include <fstream>
#include <algorithm>
#include "cmath"
#include <gtest/gtest.h>
#include <vector>

#include "Utils.hpp"
#include "DiscreteFractureNetwork.hpp"


using namespace testing;
//Test ImportFracture
// BISOGNA FARLO? E IN CASO COME LO FACCIO?



//Test ReadFracture
TEST(ReadFractureTest, FileNotFound)

{
    DiscreteFractureNetwork fracture;

    EXPECT_FALSE(ReadFracture("/nonexistent/path/", "nonexistent_file.txt", fracture));
}


TEST(ReadFractureTest, ValidFileSFracture)
{
    string filePathInput = "DFN";
    string fileName = "valid_fracture.txt";

    string content =
        "Number of Fractures: \n"
        "4\n"
        "Vertices\n"
        "X Coordinates\n"
        "0.0 1.0 2.0 3.0\n"
        "Y Coordinates\n"
        "0.0 1.0 2.0 3.0\n"
        "Z Coordinates\n"
        "0.0 1.0 2.0 3.0\n";

    {
        ofstream outFile(filePathInput + fileName);
        outFile << content;
        outFile.close();
    }

    DiscreteFractureNetwork fracture;

//    EXPECT_TRUE(ReadFracture(filePath, fileName, fracture));
//    EXPECT_EQ(fracture.numFracture, 1);
//    EXPECT_EQ(fracture.fractureID[0], 1);
//    EXPECT_EQ(fracture.NumVertices[0], 4);
//    EXPECT_EQ(fracture.vertices[0](0,0), 0.0);
//    EXPECT_EQ(fracture.vertices[0](1,1), 1.0);
//    EXPECT_EQ(fracture.vertices[0](2,2), 2.0);
//    EXPECT_EQ(fracture.vertices[0](0,3), 3.0);
}



//Test BoundingBox
TEST(BBox3DTest, BoundingBoxCalculation)
{
    MatrixXd vertices(3, 4);
    vertices << 1, 2, 3, 4,
        1, 2, 3, 4,
        1, 2, 3, 4;

    BoundingBox bbox = BBox3D(vertices);

    EXPECT_EQ(bbox.min, Vector3d(1, 1, 1));
    EXPECT_EQ(bbox.max, Vector3d(4, 4, 4));
}



// Test FractureIntersection
TEST(FractureIntersectionTest, TestSuccessfulIntersection)
{
    DiscreteFractureNetwork fracture;
    Traces trace;

    fracture.numFracture = 2;
    fracture.fractureID = {0, 1};
    fracture.vertices.resize(8);


    // Definisci i vertici della prima frattura
    fracture.vertices[0] << 0, 0, 0;
    fracture.vertices[1] << 1, 0, 0;
    fracture.vertices[2] << 1, 1, 0;
    fracture.vertices[3] << 0, 1, 0;

    // Definisci i vertici della seconda frattura
    fracture.vertices[4] << 0.8, 0, -0.1;
    fracture.vertices[5] << 0.8, 0, 0.3;
    fracture.vertices[6] << 0.8, 1, 0.3;
    fracture.vertices[7] << 0.8, 1, -0.1;


    bool intersectionResult = FractureIntersection(fracture, trace);

    EXPECT_EQ(true, intersectionResult);
}

<<<<<<< HEAD
TEST(FractureIntersectionTest, TestNoIntersection)
{
    DiscreteFractureNetwork fracture;
    Traces trace;

    fracture.numFracture = 2;
    fracture.fractureID = {0, 1};
    fracture.vertices.resize(8);

    // Definisci i vertici della prima frattura
    fracture.vertices[0] << 0, 0, 0;
    fracture.vertices[1] << 1, 0, 0;
    fracture.vertices[2] << 1, 1, 0;
    fracture.vertices[3] << 0, 1, 0;

    // Definisci i vertici della seconda frattura
    fracture.vertices[4] << -0.2, 0.5, -0.3;
    fracture.vertices[5] << 0.3, 0.5, -0.3;
    fracture.vertices[6] << 0.3, 0.5, 0.4;
    fracture.vertices[7] << -0.2, 0.5, 0.4;


    bool intersectionResult = FractureIntersection(fracture, trace);


    EXPECT_EQ(false, intersectionResult);

}


// Test SaveTraces
TEST(SaveTracesTest, CheckSize)
{
    Traces trace;
    Vector3d point(0, 0, 0);
    Vector3d s(1, 1, 1);

    double n = 1.0;
    double m = 2.0;
    unsigned int Id1 = 1;
    unsigned int Id2 = 2;

    SaveTraces(n, m, point, s, trace, Id1, Id2);

    EXPECT_EQ(trace.coordinates.size(), 1);
    EXPECT_EQ(trace.fractureId.size(), 1);
    EXPECT_EQ(trace.traceId.size(), 1);
    EXPECT_EQ(trace.length.size(), 1);
    EXPECT_EQ(trace.numTraces, 1);
}

TEST(SaveTracesTest, CheckCoordinates)
{
    Traces trace;
    Vector3d point(0, 0, 0);
    Vector3d s(1, 1, 1);

    double n = 1.0;
    double m = 2.0;
    unsigned int Id1 = 1;
    unsigned int Id2 = 2;

    SaveTraces(n, m, point, s, trace, Id1, Id2);

    EXPECT_EQ(trace.coordinates[0](0, 0), 1.0); // P.x
    EXPECT_EQ(trace.coordinates[0](1, 0), 1.0); // P.y
    EXPECT_EQ(trace.coordinates[0](2, 0), 1.0); // P.z
    EXPECT_EQ(trace.coordinates[0](0, 1), 2.0); // Q.x
    EXPECT_EQ(trace.coordinates[0](1, 1), 2.0); // Q.y
    EXPECT_EQ(trace.coordinates[0](2, 1), 2.0); // Q.z
}

//TEST(SaveTracesTest, CheckFractureAndTraceIDs)
//{
//    Traces trace;
//    Vector3d point(0, 0, 0);
//    Vector3d s(1, 1, 1);

//    double n = 1.0;
//    double m = 2.0;
//    unsigned int Id1 = 1;
//    unsigned int Id2 = 2;

//    SaveTraces(n, m, point, s, trace, Id1, Id2);

//    EXPECT_EQ(trace.fractureId, Id1);
//    EXPECT_EQ(trace.fractureId, Id2);

//    EXPECT_EQ(trace.traceId[0], Id1);
//    EXPECT_EQ(trace.traceId[0], Id2);
//}

TEST(SaveTracesTest, CheckLength)
{
    Traces trace;
    Vector3d point(0, 0, 0);
    Vector3d s(1, 1, 1);

    double n = 1.0;
    double m = 2.0;
    unsigned int Id1 = 1;
    unsigned int Id2 = 2;

    SaveTraces(n, m, point, s, trace, Id1, Id2);

    double expectedLength = PointDistance(Vector3d(1.0, 1.0, 1.0), Vector3d(2.0, 2.0, 2.0));
    EXPECT_DOUBLE_EQ(trace.length[0], expectedLength);
}



//Test TraceReorder
TEST(TraceReorderTest, NoTraces)
{
    DiscreteFractureNetwork fracture;
    Traces trace;

    fracture.numFracture = 2;
    fracture.fractureID = {0, 1};
    fracture.NumVertices = {3, 3};

    bool result = TraceReorder(fracture, trace);

    EXPECT_EQ(result, true);
    EXPECT_EQ(trace.traceReordered.size(), 2);
    EXPECT_EQ(trace.traceReordered[0].size(), 0);
    EXPECT_EQ(trace.traceReordered[1].size(), 0);
}


TEST(TraceReorderTest, TestTraceReorder)
{
    DiscreteFractureNetwork fracture;
    Traces trace;

    fracture.numFracture = 2;
    fracture.fractureID = {0, 1};
    fracture.NumVertices = {4, 4};


    fracture.vertices.resize(2);
    fracture.vertices[0].resize(3, 4);
    fracture.vertices[0] << 0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0.5, 0.5, 0;

    fracture.vertices[1].resize(3, 4);
    fracture.vertices[1] << 2, 0, 0,
        3, 0, 0,
        2, 1, 0,
        2.5, 0.5, 0;


    trace.numTraces = 4;
    trace.fractureId.resize(4);

    trace.fractureId.push_back({0});
    trace.fractureId.push_back({0});
    trace.fractureId.push_back({1});
    trace.fractureId.push_back({1});

    trace.coordinates.resize(4);
    trace.coordinates[0].resize(3, 2);
    trace.coordinates[0] << 0.1, 0.1,
        0.2, 0.2,
        0.3, 0.3;

    trace.coordinates[1].resize(3, 2);
    trace.coordinates[1] << 0.8, 0.8,
        0.9, 0.9,
        1.0, 1.0;

    trace.coordinates[2].resize(3, 2);
    trace.coordinates[2] << 2.1, 0.1,
        2.2, 0.2,
        2.3, 0.3;

    trace.coordinates[3].resize(3, 2);
    trace.coordinates[3] << 2.8, 0.8,
        2.9, 0.9,
        3.0, 1.0;

    trace.length = {0.5, 0.4, 0.6, 0.7};

    bool reorderResult = TraceReorder(fracture, trace);


    EXPECT_EQ(true, reorderResult);

    EXPECT_EQ(trace.traceReordered.size(), 2);


    vector<tuple<unsigned int, bool, double>>& reorderedTraces0 = trace.traceReordered[0];
    EXPECT_EQ(reorderedTraces0.size(), 4);


    int countPassing0 = 0;
    int countNotPassing0 = 0;
    for (size_t n = 0; n < reorderedTraces0.size(); ++n) {
        if (get<1>(reorderedTraces0[n]) == false) // Se è passante
            countPassing0++;
        else if (get<1>(reorderedTraces0[n]) == true) // Se non è passante
            countNotPassing0++;
    }
    EXPECT_EQ(countPassing0, 2);
    EXPECT_EQ(countNotPassing0, 2);

    vector<tuple<unsigned int, bool, double>>& reorderedTraces1 = trace.traceReordered[1];
    EXPECT_EQ(reorderedTraces1.size(), 0); // Verifica che non ci siano tracce riordinate per la seconda frattura
}





//Test PointDistance
//TEST(PointDistanceTest, PositiveDistance)
//{
//    Vector3d P(1.0, 2.0, 3.0);
//    Vector3d Q(4.0, 5.0, 6.0);

//    auto IsApproxEqual = [](double a, double b, double tolerance = 1e-10)
//    {
//        return abs(a - b) <= tolerance;
//    };

//    double expected_distance = sqrt(27.0);

//    EXPECT_TRUE(IsApproxEqual(PointDistance(P, Q), expected_distance));
//}


//TEST(PointDistanceTest, NegativeDistance)
//{
//    Vector3d P(1.0, 2.0, 3.0);
//    Vector3d Q(-4.0, -5.0, -6.0);

//    auto IsApproxEqual = [](double a, double b, double tolerance = 1e-6)
//    {
//        return abs(a - b) <= tolerance;
//    };

//    double expected_distance = sqrt(182.0);

//    EXPECT_TRUE(IsApproxEqual(PointDistance(P, Q), expected_distance));
//}



//Test Reordering
TEST(ReorderingTest, SameSizeVectors)
{
    vector<unsigned int> idTraces = {1, 2, 3};
    vector<double> length = {3.0, 2.0, 1.0};

    EXPECT_TRUE(reordering(idTraces, length));


    EXPECT_EQ(idTraces[0], 3);
    EXPECT_EQ(idTraces[1], 2);
    EXPECT_EQ(idTraces[2], 1);

    EXPECT_EQ(length[0], 1.0);
    EXPECT_EQ(length[1], 2.0);
    EXPECT_EQ(length[2], 3.0);\
}



//Test FindTraces
TEST(FindTracesTest, CorrectIntersection)
{
    DiscreteFractureNetwork fracture;
    Traces trace;


    fracture.numFracture = 2;
    fracture.fractureID = {0, 1};
    fracture.NumVertices = {4, 4};

    fracture.vertices.resize(2);
    fracture.vertices[0].resize(3, 4);
    fracture.vertices[0] << 0, 0, 0,
                            1, 0, 0,
                            0, 1, 0,
                            0.5, 0.5, 0;

    fracture.vertices[1].resize(3, 4);
    fracture.vertices[1] << 0.5, 0, 0,
                            1.5, 0, 0,
                            0.5, 1, 0,
                            1, 0.5, 0;



    Vector3d s(1, 1, 1);
    Vector3d point(0, 0, 0);

    bool result = FindTraces(s, point, fracture, 0, 1, trace);


    EXPECT_EQ(result, true);

    EXPECT_EQ(trace.numTraces, 1);

    EXPECT_EQ(trace.fractureId.size(), 1);
    EXPECT_EQ(trace.fractureId[0].size(), 2);

    EXPECT_EQ(trace.fractureId[0][0], 0);
    EXPECT_EQ(trace.fractureId[0][1], 1);
}


TEST(FindTracesTest, NoIntersection)
{
    DiscreteFractureNetwork fracture;
    Traces trace;


    fracture.numFracture = 2;
    fracture.fractureID = {0, 1};
    fracture.NumVertices = {4, 4};


    fracture.vertices.resize(2);
    fracture.vertices[0].resize(3, 4);
    fracture.vertices[0] << 0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0.5, 0.5, 0;

    fracture.vertices[1].resize(3, 4);
    fracture.vertices[1] << 2, 0, 0,
        3, 0, 0,
        2, 1, 0,
        2.5, 0.5, 0;


    Vector3d s(1, 1, 1);
    Vector3d point(0, 0, 0);

    bool result = FindTraces(s, point, fracture, 0, 1, trace);


    EXPECT_EQ(result, true);

    EXPECT_EQ(trace.numTraces, 0);
    EXPECT_EQ(trace.fractureId.size(), 0);
}



#endif
