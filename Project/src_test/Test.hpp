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
#include "PolygonalMesh.hpp"


using namespace testing;


//Test ReadFracture
TEST(ReadFractureTest, FileNotFound)
{
    DiscreteFractureNetwork fracture;

    EXPECT_FALSE(ReadFracture("/nonexistent/path/", "nonexistent_file.txt", fracture));
}


TEST(ReadFractureTest, ValidFileSFracture)
{
    string filePath = "DFN";
    string fileName = "/FR3_data.txt";

    DiscreteFractureNetwork fracture;

    EXPECT_TRUE(ReadFracture(filePath, fileName, fracture));

    EXPECT_EQ(fracture.numFracture, 3);
    EXPECT_EQ(fracture.fractureID[0], 0);
    EXPECT_EQ(fracture.NumVertices[0], 4);
}


//Test BoundingBox
TEST(BBox3DTest, BoundingBoxCalculation)
{
    MatrixXd vertices(3, 4);
    vertices << 0.0, 1.0, 1.0, 0.0,
                0.0, 0.0, 1.0, 1.0,
                0.0, 0.0, 0.0, 0.0;


    Vector3d expected_min(0.0, 0.0, 0.0);
    Vector3d expected_max(1.0, 1.0, 0.0);

    vector<Vector3d> bbox = BBox3D(vertices);

    EXPECT_EQ(bbox[0], expected_min);
    EXPECT_EQ(bbox[1], expected_max);
}


// Test FractureIntersection
TEST(FractureIntersectionTest, IntersectingFractures)
{
   DiscreteFractureNetwork fracture;

   fracture.numFracture = 2;
   fracture.fractureID = {0, 1};
   fracture.NumVertices = {4,4};

   MatrixXd vertices1(3, 4);
   vertices1 << 0, 1, 1, 0,
       0, 0, 1, 1,
       0, 0, 0, 0;
   MatrixXd vertices2(3, 4);
   vertices2 << -0.2378, 0.3162, 0.3162,-0.2378,
       0.5, 0.5, 0.5, 0.5,
       -0.3444, -0.3444, 0.4528, 0.4528;

   fracture.vertices = {vertices1, vertices2};

   Traces trace;
   bool result = FractureIntersection(fracture, trace);

   EXPECT_TRUE(result);
   EXPECT_EQ(trace.numTraces, 1);
   EXPECT_EQ(trace.fractureId.size(), 1);
   EXPECT_EQ(trace.fractureId[0][0] , 0);
   EXPECT_EQ(trace.fractureId[0][1] , 1);
   EXPECT_EQ(trace.coordinates.size(), 1);
}


TEST(FractureIntersectionTest, NonIntersectingFractures)
{
   DiscreteFractureNetwork fracture;
   fracture.numFracture = 2;
   fracture.fractureID = {0, 1};
   fracture.NumVertices = {4,4};


   MatrixXd vertices1(3, 4);
   vertices1 << 0.8, 0.8, 0.8, 0.8,
       0, 0, 1, 1,
       -0.1, 0.3, 0.3, -0.1;

   MatrixXd vertices2(3, 4);
   vertices2 << -0.2, 0.3, 0.3, -0.2,
                 0.5, 0.5, 0.5, 0.5,
                -0.3, -0.3, 0.4, 0.4;

   fracture.vertices = {vertices1, vertices2};

   Traces trace;
   bool result = FractureIntersection(fracture, trace);

   EXPECT_TRUE(result);
   EXPECT_EQ(trace.numTraces, 0);
}


//Test FindTraces
TEST(FindTracesTest, Intersection)
{
    DiscreteFractureNetwork fracture;
    Traces trace;

    fracture.numFracture = 2;
    fracture.fractureID = {0, 1};
    fracture.NumVertices = {4, 4};

    fracture.vertices.resize(2);


    fracture.vertices[0].resize(3, 4);
    fracture.vertices[0] << 0, 1, 1, 0,
        0, 0, 1, 1,
        0, 0, 0, 0;

    fracture.vertices[1].resize(3, 4);
    fracture.vertices[1] << -0.2378, 0.3162, 0.3162,-0.2378,
        0.5, 0.5, 0.5, 0.5,
        -0.3444, -0.3444, 0.4528, 0.4528;


    Vector3d s(0.4035, 0, 0);
    Vector3d point(0, 0.5, 0);

    bool result = FindTraces(s, point, fracture, 0, 1, trace);

    EXPECT_TRUE(result);
    EXPECT_EQ(trace.numTraces, 1);
    EXPECT_EQ(trace.fractureId.size(), 1);
    EXPECT_EQ(trace.coordinates.size(), 1);
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
    fracture.vertices[0] << 0, 1, 0.5, 0,
                            0, 0, 0.5, 1,
                            0, 0, 0, 0;

    fracture.vertices[1].resize(3, 4);
    fracture.vertices[1] << 2, 3, 2.5, 2,
                            0, 0, 0.5, 1,
                            0, 0, 0,   0;


    Vector3d s(1, 1, 1);
    Vector3d point(0, 0, 0);

    bool result = FindTraces(s, point, fracture, 0, 1, trace);


    EXPECT_EQ(result, true);

    EXPECT_EQ(trace.numTraces, 0);
    EXPECT_EQ(trace.fractureId.size(), 0);
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

    EXPECT_EQ(trace.coordinates[0](0, 0), 2.0); // P.x
    EXPECT_EQ(trace.coordinates[0](1, 0), 2.0); // P.y
    EXPECT_EQ(trace.coordinates[0](2, 0), 2.0); // P.z
    EXPECT_EQ(trace.coordinates[0](0, 1), 1.0); // Q.x
    EXPECT_EQ(trace.coordinates[0](1, 1), 1.0); // Q.y
    EXPECT_EQ(trace.coordinates[0](2, 1), 1.0); // Q.z
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

    Vector3d P = point + s * m;
    Vector3d Q = point + s * n;

    double expectedLength = sqrt(((P[0]-Q[0])*(P[0]-Q[0])+((P[1]-Q[1])*(P[1]-Q[1]))+((P[2]-Q[2])*(P[2]-Q[2]))));
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
        1, 1, 0,
        0, 1, 0;

    fracture.vertices[1].resize(3, 4);
    fracture.vertices[1] << 2, 0, 0,
        3, 0, 0,
        3, 1, 0,
        2, 1, 0;

    trace.numTraces = 4;
    trace.fractureId = {{0, 0}, {0, 0}, {1, 1}, {1, 1}};
    trace.traceId = {0, 1, 2, 3};

    trace.coordinates.resize(4);
    trace.coordinates[0].resize(2, 3);
    trace.coordinates[0] << 0.5, 0.1, 0,
        0.5, 0.9, 0;

    trace.coordinates[1].resize(2, 3);
    trace.coordinates[1] << 0.8, 0.1, 0,
        0.9, 0.8, 0;

    trace.coordinates[2].resize(2, 3);
    trace.coordinates[2] << 2.5, 0.1, 0,
        2.5, 0.9, 0;

    trace.coordinates[3].resize(2, 3);
    trace.coordinates[3] << 2.8, 0.1, 0,
        2.9, 0.8, 0;

    trace.length = {0.5, 0.4, 0.6, 0.7};

    bool reorderResult = TraceReorder(fracture, trace);

    EXPECT_TRUE(reorderResult);
    EXPECT_EQ(trace.traceReordered.size(), 2);

    for (size_t i = 0; i < trace.traceReordered.size(); ++i) {
        vector<tuple<unsigned int, bool, double>>& reorderedTraces = trace.traceReordered[i];
        EXPECT_EQ(reorderedTraces.size(), 2);


        int countPassing = 0;
        int countNotPassing = 0;
        for (const auto& reorderedTrace : reorderedTraces)
        {
            if (get<1>(reorderedTrace) == false)
                countPassing++;

            else if (get<1>(reorderedTrace) == true)
                countNotPassing++;
        }
        EXPECT_EQ(countPassing, 0);
        EXPECT_EQ(countNotPassing, 2);
    }
}



//Test Reordering
TEST(ReorderingTest, SameSizeVectors)
{
    vector<unsigned int> idTraces = {1, 2, 3};
    vector<double> length = {3.0, 2.0, 1.0};

    EXPECT_TRUE(reordering(idTraces, length));


    EXPECT_EQ(idTraces[0], 1);
    EXPECT_EQ(idTraces[1], 2);
    EXPECT_EQ(idTraces[2], 3);

    EXPECT_EQ(length[0], 3.0);
    EXPECT_EQ(length[1], 2.0);
    EXPECT_EQ(length[2], 1.0);
}
// *****************************************************************************************************************************************************



//Test createSubfracture
TEST(CreateSubfractureTest, CreationSubfracture)
{
    vector<Vector3d> subfracture = {
        Vector3d(0.0, 0.0, 0.0),
        Vector3d(1.0, 0.0, 0.0),
        Vector3d(1.0, 1.0, 0.0),
        Vector3d(0.0, 1.0, 0.0)
    };


    vector<Vector3d> cuttingTrace = {
        Vector3d(0.5, 0.0, 0.0),
        Vector3d(0.5, 1.0, 0.0)
    };


    vector<vector<Vector3d>> subfractureVertices1;


    bool result = createSubfracture(subfracture, cuttingTrace, subfractureVertices1);


    ASSERT_TRUE(result);


    ASSERT_EQ(subfractureVertices1.size(), 2);


    vector<Vector3d> expectedSubfracture1 = {
        Vector3d(0.0, 0.0, 0.0),
        Vector3d(0.5, 0.0, 0.0),
        Vector3d(0.5, 1.0, 0.0),
        Vector3d(0.0, 1.0, 0.0)
    };

    ASSERT_EQ(subfractureVertices1[0].size(), expectedSubfracture1.size());

    vector<Vector3d> expectedSubfracture2 = {
        Vector3d(0.5, 0.0, 0.0),
        Vector3d(1.0, 0.0, 0.0),
        Vector3d(1.0, 1.0, 0.0),
        Vector3d(0.5, 1.0, 0.0)
    };

    ASSERT_EQ(subfractureVertices1[1].size(), expectedSubfracture2.size());
}



//Test extendTraces
TEST(ExtendTracesTest, CorrectIntersectionPoints)
{
    vector<Vector3d> subFractureVertices = {
        Vector3d(0, 0, 0),
        Vector3d(1, 0, 0),
        Vector3d(1, 1, 0),
        Vector3d(0, 1, 0)
    };

    vector<Vector3d> subTraceVertices = {
        Vector3d(0.5, 0, 0),
        Vector3d(0.5, 1, 0)
    };

    vector<Vector3d> expectedExtendedVertices = {
        Vector3d(0.5, 0, 0),
        Vector3d(0.5, 1, 0)
    };

    vector<Vector3d> result = extendTraces(subFractureVertices, subTraceVertices);

    EXPECT_EQ(result.size(), expectedExtendedVertices.size());
    for (size_t i = 0; i < result.size(); ++i) {
        EXPECT_EQ(result[i].x(), expectedExtendedVertices[i].x());
        EXPECT_EQ(result[i].y(), expectedExtendedVertices[i].y());
        EXPECT_EQ(result[i].z(), expectedExtendedVertices[i].z());
    }
}



//Test createMesh
TEST(CreateMeshTest, CorrectMeshCreation)
{
    vector<pair<vector<Vector3d>, vector<pair<vector<Vector3d>, unsigned int>>>> subFracture = {
        {
            {
                Vector3d(0, 0, 0),
                Vector3d(1, 0, 0),
                Vector3d(1, 1, 0),
                Vector3d(0, 1, 0)
            },
            {}
        }
    };

    PolygonalMesh mesh;

    bool result = createMesh(subFracture, mesh);

    ASSERT_TRUE(result);

    // Verifica che le coordinate dei vertici siano corrette
    ASSERT_EQ(mesh.coordinates0D.size(), 4);
    EXPECT_EQ(mesh.coordinates0D[0], Vector3d(0, 0, 0));
    EXPECT_EQ(mesh.coordinates0D[1], Vector3d(1, 0, 0));
    EXPECT_EQ(mesh.coordinates0D[2], Vector3d(1, 1, 0));
    EXPECT_EQ(mesh.coordinates0D[3], Vector3d(0, 1, 0));

    // Verifica che gli ID delle celle 0D siano corretti
    ASSERT_EQ(mesh.cellId0D.size(), 4);
    EXPECT_EQ(mesh.cellId0D[0], 0);
    EXPECT_EQ(mesh.cellId0D[1], 1);
    EXPECT_EQ(mesh.cellId0D[2], 2);
    EXPECT_EQ(mesh.cellId0D[3], 3);

    // Verifica che gli ID delle celle 1D siano corretti
    ASSERT_EQ(mesh.cellId1D.size(), 4);
    EXPECT_EQ(mesh.cellId1D[0], 0);
    EXPECT_EQ(mesh.cellId1D[1], 1);
    EXPECT_EQ(mesh.cellId1D[2], 2);
    EXPECT_EQ(mesh.cellId1D[3], 3);

    // Verifica che gli ID dei vertici delle celle 1D siano corretti
    ASSERT_EQ(mesh.verticesId1D.size(), 4);
    EXPECT_EQ(mesh.verticesId1D[0], Vector2i(0, 1));
    EXPECT_EQ(mesh.verticesId1D[1], Vector2i(1, 2));
    EXPECT_EQ(mesh.verticesId1D[2], Vector2i(2, 3));
    EXPECT_EQ(mesh.verticesId1D[3], Vector2i(3, 0));

    // Verifica che gli ID delle celle 2D siano corretti
    ASSERT_EQ(mesh.cellId2D.size(), 1);
    EXPECT_EQ(mesh.cellId2D[0], 0);

    // Verifica che il numero di vertici delle celle 2D sia corretto
    ASSERT_EQ(mesh.numVertices2D.size(), 1);
    EXPECT_EQ(mesh.numVertices2D[0], 4);

    // Verifica che gli ID dei vertici delle celle 2D siano corretti
    ASSERT_EQ(mesh.verticesId2D.size(), 1);
    ASSERT_EQ(mesh.verticesId2D[0].size(), 4);
    EXPECT_EQ(mesh.verticesId2D[0][0], 0);
    EXPECT_EQ(mesh.verticesId2D[0][1], 1);
    EXPECT_EQ(mesh.verticesId2D[0][2], 2);
    EXPECT_EQ(mesh.verticesId2D[0][3], 3);

    // Verifica che il numero di lati delle celle 2D sia corretto
    ASSERT_EQ(mesh.numEdges2D.size(), 1);
    EXPECT_EQ(mesh.numEdges2D[0], 4);

    // Verifica che gli ID dei lati delle celle 2D siano corretti
    ASSERT_EQ(mesh.edgesId2D.size(), 1);
    ASSERT_EQ(mesh.edgesId2D[0].size(), 4);
    EXPECT_EQ(mesh.edgesId2D[0][0], 0);
    EXPECT_EQ(mesh.edgesId2D[0][1], 1);
    EXPECT_EQ(mesh.edgesId2D[0][2], 2);
    EXPECT_EQ(mesh.edgesId2D[0][3], 3);
}



#endif
