/**
 * Project: Interactive ARAP
 * File:    TestArap.h
 * Authors: Ankur Deria, Michael Dey, Bendeguz Timar, Andrea Solanas de Vicente
 */

#pragma once

#include "Mesh.h"
#include <igl/readOFF.h>
#include <igl/upsample.h>
#include <chrono>
#include "ArapAlgorithm.h"

 /*
  * Class that implements tests for the arap algorithm class.
  */
struct TestArap {
public:
	// constructor 1
	TestArap(std::map<int, bool> _fixedFaces);
	// Constrcutor 2
	TestArap(const std::string& path);

	std::shared_ptr<Eigen::MatrixXd> vertices; // Vertices from the mesh
	std::shared_ptr<Eigen::MatrixXi> faces; // Faces from the mesh
	int moving_vertex = -1; // moving vertex
	std::map<int, bool> fixed_faces_selected; // map of <face_id, true/false> inc ase it is selected
	ArapAlgorithm aa = ArapAlgorithm(); // arap algorithm instance for testing

	int testIndex = 0;

private:
	void ExecuteTest(); // Function to execute test
};
