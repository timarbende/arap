/**
 * Project: Interactive ARAP
 * File:    main.cpp
 * Authors: Ankur Deria, Michael Dey, Bendeg�z Timar, Andrea Solanas de Vicente
 */

#include <ArapAlgorithm.h>
#include <igl/readOFF.h>;
#include <igl/opengl/glfw/Viewer.h>;
#include <igl/unproject_onto_mesh.h>;
#include <Mesh.h>;
#include <TestArap.h>;

const bool test = false;

/*
 * Main function
 */
int main(int argc, char** argv)
{
    std::string filePath = "../meshes/dino/dino.off"; // filepath where mesh can be found
    //std::string filePath = "../meshes/test/octahedron_test.off";
    if (!test) {
        Mesh mesh(filePath); // create the mesh
        mesh.launchViewer(); // launch the viewer
    }
    else {
        // Choose the fixed faces, can be random
        TestArap test1 = TestArap(filePath);
    }

    return EXIT_SUCCESS;
}

