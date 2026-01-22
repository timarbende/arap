/**
 * Project: Interactive ARAP
 * File:    TestArap.cpp
 * Authors: Ankur Deria, Michael Dey, Bendegúz Timár, Andrea Solanas de Vicente
 */

#include "TestArap.h"

 /**
  * Test constructor 2
  * 
  * Executes the test with the example stated inside the function. (Octahedron)
  *
  * @param _fixedFaces: map of the ids of the fixed faces and if they are fixed or not
  * @return TestArap instance
  */
TestArap::TestArap(std::map<int, bool> _fixedFaces) : fixed_faces_selected(_fixedFaces) {

    static int testIndex = 0;
    testIndex++;
    std::cout << "===== Test n " << testIndex << ": =====" << std::endl;

    // Octahedron
    vertices = std::make_shared<Eigen::MatrixXd>((Eigen::MatrixXd(6, 3) <<
        0.0, 0.0, 2.0,
        2.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        - 2.0, 0.0, 0.0,
        0.0, - 2.0, 0.0,
        0.0, 0.0, - 2.0).finished());
    faces = std::make_shared<Eigen::MatrixXi>((Eigen::MatrixXi(8, 3) <<
        1, 0, 4,
        4, 0, 3,
        3, 0, 2,
        2, 0, 1,
        1, 5, 2,
        2, 5, 3,
        3, 5, 4,
        4, 5, 1).finished().array());

    // Complete the meshes initialization & Start the Test
    ExecuteTest();
}

/**
  * Test constructor 1
  *
  * Executes the test with the example passed as file path
  *
  * @param _fixedFaces: map of the ids of the fixed faces and if they are fixed or not
  * @param path: filepath where the off file for the mesh is stored
  * @return TestArap instance
  */
TestArap::TestArap(const std::string& path) {
    testIndex++;
    std::cout << "===== Test n " << testIndex << ": =====" << std::endl;

    // Load a mesh in OFF format
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    if (!igl::readOFF(path, V, F)){
        DEBUG_STDERR("Failed to read file: ");
        DEBUG_STDERR(path);
    }
    vertices = std::make_shared<Eigen::MatrixXd>(std::move(V));
    faces = std::make_shared<Eigen::MatrixXi>(std::move(F));

    // Complete the meshes initialization & Start the Test
    ExecuteTest();
}

/**
  * Execute Test
  *
  * Executes the arap algorithm by manually setting the fixed faces and moving vertex. Checks if data is correct
  *
  * @return void
  */
void TestArap::ExecuteTest() {
    // random fixed faces
    std::vector<int> fixed_Faces = { 4,5 };

    //choose moving vertex
    moving_vertex = 0;

    // precompute arap
    aa.preComputeDeformation(vertices, faces);

    // compute deformation
    Eigen::Vector3f position_vertex = Eigen::Vector3f(0.0, 0.0, 2.0); // change the position if wanted
    aa.updateMovingVertex(moving_vertex, position_vertex, faces, fixed_Faces);

    // call arap algorithm
    std::shared_ptr<Eigen::MatrixXd> newVertices = aa.doArap(vertices); 

    std::cout << "===== Test n " << testIndex << ": FINISHED ======" << std::endl;
}