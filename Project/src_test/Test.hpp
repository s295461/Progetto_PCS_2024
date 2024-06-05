#ifndef TEST_H
#define TEST_H

#include <iostream>
#include "Eigen/Eigen"
#include <fstream>
#include <algorithm>
#include "cmath"
#include <gtest/gtest.h>

#include "Utils.hpp"
#include "DiscreteFractureNetwork.hpp"


using namespace testing;
using namespace FractureNetwork;


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


//Test PointDistance
TEST(PointDistanceTest, PositiveDistance)
{
    Vector3d P(1.0, 2.0, 3.0);
    Vector3d Q(4.0, 5.0, 6.0);

    auto IsApproxEqual = [](double a, double b, double tolerance = 1e-10)
    {
        return abs(a - b) <= tolerance;
    };

    double expected_distance = sqrt(27.0);

    EXPECT_TRUE(IsApproxEqual(PointDistance(P, Q), expected_distance));
}


TEST(PointDistanceTest, NegativeDistance)
{
    Vector3d P(1.0, 2.0, 3.0);
    Vector3d Q(-4.0, -5.0, -6.0);

    auto IsApproxEqual = [](double a, double b, double tolerance = 1e-6)
    {
        return abs(a - b) <= tolerance;
    };

    double expected_distance = sqrt(182.0);

    EXPECT_TRUE(IsApproxEqual(PointDistance(P, Q), expected_distance));
}



//Test Reordering
TEST(ReorderingTest, SameSizeVectors) {
    vector<unsigned int> idTraces = {1, 2, 3};
    vector<double> length = {3.0, 2.0, 1.0};

    EXPECT_TRUE(reordering(idTraces, length));


    EXPECT_EQ(idTraces[0], 3);
    EXPECT_EQ(idTraces[1], 2);
    EXPECT_EQ(idTraces[2], 1);

    EXPECT_EQ(length[0], 1.0);
    EXPECT_EQ(length[1], 2.0);
    EXPECT_EQ(length[2], 3.0);
}



//Test ReadFracture
TEST(ReadFractureTest, FileNotFound)
{
    DiscreteFractureNetwork fracture;

    EXPECT_FALSE(ReadFracture("/nonexistent/path/", "nonexistent_file.txt", fracture));
}


TEST(ReadFractureTest, ValidFileSingleFracture)
{
    string filePath = "./";
    string fileName = "valid_single_fracture.txt";

    string content =
        "Number of Fractures: 1\n"
        "1 4\n"
        "Vertices\n"
        "X Coordinates\n"
        "0.0 1.0 2.0 3.0\n"
        "Y Coordinates\n"
        "0.0 1.0 2.0 3.0\n"
        "Z Coordinates\n"
        "0.0 1.0 2.0 3.0\n";

    {
        ofstream outFile(filePath + fileName);
        outFile << content;
        outFile.close();
    }

    DiscreteFractureNetwork fracture;

    EXPECT_TRUE(ReadFracture(filePath, fileName, fracture));
    EXPECT_EQ(fracture.numFracture, 1);
    EXPECT_EQ(fracture.fractureID[0], 1);
    EXPECT_EQ(fracture.NumVertices[0], 4);
    EXPECT_EQ(fracture.vertices[0](0,0), 0.0);
    EXPECT_EQ(fracture.vertices[0](1,1), 1.0);
    EXPECT_EQ(fracture.vertices[0](2,2), 2.0);
    EXPECT_EQ(fracture.vertices[0](0,3), 3.0);
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

TEST(SaveTracesTest, CheckFractureAndTraceIDs)
{
    Traces trace;
    Vector3d point(0, 0, 0);
    Vector3d s(1, 1, 1);

    double n = 1.0;
    double m = 2.0;
    unsigned int Id1 = 1;
    unsigned int Id2 = 2;

    SaveTraces(n, m, point, s, trace, Id1, Id2);

    EXPECT_EQ(trace.fractureId , Id1);
    EXPECT_EQ(trace.fractureId , Id2);
    EXPECT_EQ(trace.traceId[0], 0);
}

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

#endif
