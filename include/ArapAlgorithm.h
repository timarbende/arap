/**
 * Project: Interactive ARAP
 * File:    ArapAlgorithm.h
 * Authors: Ankur Deria, Michael Dey, Bendeguz Timar, Andrea Solanas de Vicente
 */
#pragma once

#include <Eigen/Dense>;
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/StdVector>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/NonLinearOptimization>
#include <cstdlib>;
#include <iostream>
#include <map>
#include <chrono>
#ifdef OMP
#include <omp.h>
#endif

#ifdef _DEBUG //Use this method to print all statements and throw errors so that it only works in debug builds and not in release builds
#define DEBUG_STDERR(x) (std::cerr << (x) << std::endl)
#define DEBUG_STDOUT(x) (std::cout << (x) << std::endl)
#else 
#define DEBUG_STDERR(x)
#define DEBUG_STDOUT(x)
#endif

/* 
 * Class that implements the ARAP algorithm.
 */  
class ArapAlgorithm {
public:
	// Constructor for the class
	ArapAlgorithm(); 
	// Destructor for the class
	~ArapAlgorithm() = default;

	// Function to computer the ARAP algorithm
	std::shared_ptr<Eigen::MatrixXd> doArap(std::shared_ptr<Eigen::MatrixXd> vertices); 
	// Function to precompute the deformation of the mesh
	void preComputeDeformation(std::shared_ptr<Eigen::MatrixXd> vertices, std::shared_ptr<Eigen::MatrixXi> faces);
	// Function to update the moving vertex
	void updateMovingVertex(int movingVertex, const Eigen::Vector3f& movingVertexPosition, std::shared_ptr<Eigen::MatrixXi> faces, const std::vector<int>& fixedFaces);

private:

	// VARIABLES
	std::shared_ptr<Eigen::MatrixXd> verticesPreDeformation; // Keep the vertices pre deformation fixed for solving the linear system
	std::map<int, std::vector<int>> neighbors; // Neighborhood of vertices (Mapping between vertex id and its neighbor ids)
	Eigen::MatrixXd weight_Matrix; 	// Weight matrix used
	Eigen::MatrixXd system_Matrix_L; // System matrix -> L
	int moving_Vertex{}; // ARAP variables (the moving vertex and its position)
	Eigen::Vector3d movingVertex_Position{}; // Moving vertex position
	std::vector<int> fixed_Vertices; // -> chosen by user, the fixed vertices
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver; // USE THIS TO SOLVE and compute the deformation

	// FUNCTIONS
	// Function to obtain the neighbors of the vertices
	void obtainNeighbors(std::shared_ptr<Eigen::MatrixXi> faces);
	// Function to initialize the weights matrix
	void InitWeightMatrix(std::shared_ptr<Eigen::MatrixXi> faces);
	// Function to compute the system matrix, also called L
	void compute_L();
	// Function to obtain the fixed vertices
	void obtainFixedVertices(std::shared_ptr<Eigen::MatrixXi> faces, const std::vector<int>& fixedFaces);
	// Function to update the system matrix on fixed vertices
	void update_L_After_FixedVertices();
	// Function used to calculate the right hand side of the equation -> b
	// system L * p' = b -> only calculate b
	std::shared_ptr <Eigen::MatrixXd> compute_b(std::vector<Eigen::Matrix3d> rotationMatrices);
	// Function to calculate the rotation matrix
	// here the covariance matrix is also calculated
	std::vector<Eigen::Matrix3d> calculateRotationMatrix(std::shared_ptr<Eigen::MatrixXd> changedVertices);
	// Function to calculate the energy
	double calculateEnergy(std::shared_ptr<Eigen::MatrixXd> changedVertices, std::vector<Eigen::Matrix3d> rotationMatrices);
};
