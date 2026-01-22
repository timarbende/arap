/**
 * Project: Interactive ARAP
 * File:    ArapAlgorithm.cpp
 * Authors: Ankur Deria, Michael Dey, Bendegúz Timár, Andrea Solanas de Vicente
 */
#include <ArapAlgorithm.h>

 /**
   * Arap algorithm Constructor
   *
   * Constructs the arap algorithm instance
   *
   * @return ArapAlgorithm instance
   */
ArapAlgorithm::ArapAlgorithm() {
#ifdef OMP  // if we have OpenMP installed
    // Init Eigen in parallel mode
    Eigen::initParallel();

    // Determine the number of threads
    int nThreads;
#pragma omp parallel default(none) shared(nThreads)
    nThreads = omp_get_thread_num(); // use OpenMP to determine number of threads

    // Set number of threads Eigen
	// (Eigen uses almost 100% of CPU capacity)
    Eigen::setNbThreads(nThreads / 2);
#endif
}

/**
  * Obtain the neighbors
  *
  * Calculates the neighbors of the vertices in the mesh.
  *
  * @param faces: faces of the mesh
  * @return void
  */
void ArapAlgorithm::obtainNeighbors(std::shared_ptr<Eigen::MatrixXi> faces) {
	// Clear before starting
	neighbors.clear();
	// obtain the data since they are private values
	auto& new_neighbors = neighbors;
	auto _initial_Vertices = verticesPreDeformation;
	auto& _faces = *faces;
	
	// use OpenMP
#ifdef OMP
#pragma omp parallel for default(none) \
            shared(_initial_Vertices, new_neighbors, faces,_faces)
#endif
	// iterate over all the vertices
	for (int vertex = 0; vertex < _initial_Vertices->rows(); vertex++) {
		// vectors to save the specific neighbors 
		std::vector<int> allNeighbors;
		// Iterate over the faces to find all the neighbors
		for (int face = 0; face < faces->rows(); face++) {
			for (int j = 0; j < 3; j++) { // Iterate over the vertices (3)
				int test = _faces(face, j);
				if (test == vertex) {
					// If the face has the vertex id of this vertex (neighbor)
					test = _faces(face, (j + 1) % 3);
					allNeighbors.push_back(test); // add the next column in face matrix (neighbor 1) 
					test = _faces(face, (j + 2) % 3);
					allNeighbors.push_back(test); // add the previous column in face matrix (neighbor 2)
				}
			}

		}

		// order the vector from 0 to n
		sort(allNeighbors.begin(), allNeighbors.end());
		// eliminate those repeated vertices
		allNeighbors.erase(unique(allNeighbors.begin(), allNeighbors.end()), allNeighbors.end());

#pragma omp critical
		// Add the neighbors to the vertex -> this area is critical in parallel
		new_neighbors[vertex] = allNeighbors;
	}
}

/**
  * Initialize the weight matrix
  *
  * Initializes the weight matrix using the cotangent weights formula
  *
  * @param faces: faces of the mesh
  * @return void
  */
void ArapAlgorithm::InitWeightMatrix(std::shared_ptr<Eigen::MatrixXi> faces)
{
	// Contangent weights from here on
	weight_Matrix = Eigen::MatrixXd::Zero(verticesPreDeformation->rows(), verticesPreDeformation->rows());
	std::vector<std::vector<Eigen::Vector2i>> vertexEdges(verticesPreDeformation->rows());

	auto undeformedVertices = verticesPreDeformation;
	auto& neighborhood = neighbors;
	auto& weightMatrix = weight_Matrix;

#ifdef OMP
#pragma omp parallel for default(none) \
            shared(undeformedVertices, faces, vertexEdges)
#endif
	
	// Get the edges connected to each vertex 
	for (int i = 0; i < undeformedVertices->rows(); i++) { // Iterate over the vertices
		std::vector<Eigen::Vector2i> edges{};
		for (int j = 0; j < faces->rows(); j++) {	// Iterate over the faces
			Eigen::Vector3i face = faces->row(j);
			for (int k = 0; k < 3; k++)	{	// Iterate over the triangle
				if (face[k] == i) {
					edges.emplace_back(face[(k + 1) % 3], face[(k + 2) % 3]);
				}
			}
		}
		vertexEdges[i] = edges; // store the edge inside the vertex edges
	}

#ifdef OMP
#pragma omp parallel for default(none) \
            shared(undeformedVertices, neighborhood, weightMatrix, vertexEdges)
#endif

	for (int i = 0; i < undeformedVertices->rows(); i++) { // Iterate over the vertices
		for (int neighbor : neighborhood[i]) { // Iterate over the neighbors
			double totalAngle = 0.0;
			for (const Eigen::Vector2i& edge : vertexEdges[i]) { // Iterate over the edges
				double normBC = (undeformedVertices->row(edge[0]) - undeformedVertices->row(edge[1])).norm(); // Norm between B and C
				double normAC = (undeformedVertices->row(i) - undeformedVertices->row(edge[1])).norm(); // Norm between A and C
				double normAB = (undeformedVertices->row(i) - undeformedVertices->row(edge[0])).norm(); // Norm between A and B

				// From cosine law
				double beta = acos(((normAB * normAB) + (normBC * normBC) - (normAC * normAC)) / (2 * normAB * normBC));

				// Add to total angle if one of the points on the edge is the current neighbor
				totalAngle += (edge[0] == neighbor) * abs(tan(M_PI_2 - beta));
				totalAngle += (edge[1] == neighbor) * abs(tan(M_PI_2 - beta));
			}
			weightMatrix(i, neighbor) = abs(totalAngle) / 2;
		}
		weightMatrix(i, i) = 1.0; // Override the diagonal entry
	}
}

/**
  * Compute system matrix
  *
  * Compute the system matric using the formula obtained from paper
  * SUM over all neighbors: wij*(p'i - p'j)
  *
  * @return void
  */
void ArapAlgorithm::compute_L() {
	// initialize the system matrix L
	// size L[pre_deformed_vertices, pre_deformed_vertices]
	system_Matrix_L = Eigen::MatrixXd::Zero(verticesPreDeformation->rows(), verticesPreDeformation->rows());

	auto undeformedVertices = verticesPreDeformation;
	auto& neighborhood = neighbors;
	auto& weightMatrix = weight_Matrix;
	auto& systemMatrix = system_Matrix_L;

#ifdef OMP
#pragma omp parallel for default(none) \
            shared(undeformedVertices, neighborhood, weightMatrix, systemMatrix)
#endif

	// Iterate over the vertices
	for (int i = 0; i < undeformedVertices->rows(); i++) { 
		// iterate over the neighbors of this vertex
		for (int j : neighborhood[i]) { 
			// Follow formula in paper
			// SUM over all neighbors: wij*(p'i - p'j)
			// in the weight matrix we have the weight between i and j 
			systemMatrix(i, i) += weightMatrix(i, j); // in the diagonal
			// in the rest 
			systemMatrix(i, j) -= weightMatrix(i, j);
		}
	}
}

/**
  * Pre compute the deformation
  *
  * Compute the precomputation of the deformation. This includes obtaining the neighbors of the mesh, initializing the weight 
  * matrix and computing the initial system matrix
  *
  * @param vertices: vertices of the mesh
  * @param faces: faces of the mesh
  * @return void
  */
void ArapAlgorithm::preComputeDeformation(std::shared_ptr<Eigen::MatrixXd> vertices, std::shared_ptr<Eigen::MatrixXi> faces) {
	// Copy the vertices to the initial matrix before any deformation operation
	verticesPreDeformation = std::make_shared<Eigen::MatrixXd>(*vertices);
	// obtain the neighbors
	obtainNeighbors(faces);
	// weights should go here
	InitWeightMatrix(faces);
	// system matrix -> L (left side of system)
	compute_L(); // L p' = b -> computes only L
}

/**
  * Obtain the fixed vertices
  *
  * Calculates the fixed vertices through the fixed faces chosen by the user.
  *
  * @param faces: faces of the mesh
  * @param fixedFaces: fixed faces chosen by user 
  * @return void
  */
void ArapAlgorithm::obtainFixedVertices(std::shared_ptr <Eigen::MatrixXi> faces, const std::vector<int>& fixedFaces) {
	// Clear everything before use
	fixed_Vertices.clear();
	// save memory for the fixed vertices
	fixed_Vertices.reserve(fixedFaces.size() * faces->row(0).cols() + 1);
	// The vertices of a specific face
	Eigen::VectorXi faceVertices;
	// Add the vertices of each face to the fixed vertices
	for (int ff : fixedFaces) {
		// save the fixed face from the matrix faces into the faceVertices
		faceVertices = faces->row(ff);
		// go over the columns in faces (should be 3)
		for (int i = 0; i < 3; i++) {
			// add the vertices from faces 
			fixed_Vertices.push_back(faceVertices(i));
		}
	}
	// Add the selected vertex to the fixed vertices
	fixed_Vertices.push_back(moving_Vertex);
}

/**
  * Update system matrix after fixed vertices
  *
  * Updates the system matrix again after the fixed vertices have been chosen by the user
  *
  * @return void
  */
void ArapAlgorithm::update_L_After_FixedVertices() {
	// Update system matrix on fixed vertices to keep the fixed vertices stay where they are
	for (int fixedVertex : fixed_Vertices) { // Iterate over the fixed vertices
		//std::cout << "Number of vertex that is fixed: " << fixedVertex << std::endl;
		// in the system matrix set to zero the whole row
		system_Matrix_L.row(fixedVertex).setZero();
		// and set the diagonal for the multiplication to 1
		system_Matrix_L(fixedVertex, fixedVertex) = 1;
	}
}

/**
  * Update the moving vertex
  *
  * Updates the moving vertex to the one selected by the user. Also calls for obtaining the fixed vertices 
  * and updating the system matrix
  *
  * @param movingVertex: id of the moving vertex
  * @param movingVertexPosition: 3d position of the moving vertex
  * @param faces: faces of the mesh
  * @param fixedFaces: fixed faces selected by user
  * @return void
  */
void ArapAlgorithm::updateMovingVertex(int movingVertex, const Eigen::Vector3f& movingVertexPosition, std::shared_ptr <Eigen::MatrixXi> faces, const std::vector<int>& fixedFaces) {
	//Save the new moving vertex and its position 
	moving_Vertex = movingVertex;
	movingVertex_Position = movingVertexPosition.cast<double>();
	// now we obtain the fixed vertices (from user?)
	obtainFixedVertices(faces, fixedFaces);
	// once we have the vertices we update the sytem matrix L
	update_L_After_FixedVertices();
}

/**
  * Calculate the rotation matrices
  *
  * Compute the rotation matrices of all the vertices except the ones not moving
  *
  * @param changedVertices: matrix with the deformed matrices
  * @return vector of Matrix3D with the rotation matrix for every vertex
  */
std::vector<Eigen::Matrix3d> ArapAlgorithm::calculateRotationMatrix(std::shared_ptr<Eigen::MatrixXd> changedVertices) {
	// Create a vector to store all the rotation matrixes for each vertex
	std::vector<Eigen::Matrix3d> all_RotationMatrices(verticesPreDeformation->rows());
	// Save data since they are private values
	auto _initial_Vertices = verticesPreDeformation;
	auto& allNeighbors = neighbors;
	auto& weightMatrix = weight_Matrix;

#ifdef OMP
#pragma omp parallel for default(none) \
            shared(_initial_Vertices, allNeighbors, weightMatrix, changedVertices, all_RotationMatrices)
#endif

	for (int i = 0; i < _initial_Vertices->rows(); i++) {
		// save the number of neightbors per vertex
		const long numNeighbors = (long)(allNeighbors[i].size());
		// The definitions for the matrices P, D and P_prime can be found in the paper!
		// P: Pi is the 3 x |N(vi)| containing eij as its columns (eij : = pi − pj)
		// D: diagonal matrix containing the weights wij
		// P_prime: similar to P
		Eigen::MatrixXd P = Eigen::MatrixXd::Zero(3, numNeighbors);
		Eigen::MatrixXd D = Eigen::MatrixXd::Zero(numNeighbors, numNeighbors);
		Eigen::MatrixXd P_prime = Eigen::MatrixXd::Zero(3, numNeighbors);
		// fill the matrixes
		for (int j = 0; j < numNeighbors; j++) { // Iterate over the neighbors
			// for every column in P : take the vertex in initial vertices and substract the vertex that is neighbor of this vertex
			// pi - pj
			P.col(j) = _initial_Vertices->row(i) - _initial_Vertices->row(allNeighbors[i][j]);
			// the diagonal of D is the weight ij
			D(j, j) = weightMatrix(i, allNeighbors[i][j]);
			// Take the changed vertices now
			P_prime.col(j) = changedVertices->row(i) - changedVertices->row(allNeighbors[i][j]);
		}
		// S, the covariance matrix ( from paper)
		Eigen::Matrix3d S = P * D * P_prime.transpose();
		// using this we obtain the SVD of the covariance matrix to get the rotation matrix (from paper)
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
		const Eigen::Matrix3d& U = svd.matrixU();
		const Eigen::Matrix3d& V = svd.matrixV();

		// Computation of matrix I is necessary since UV' is only orthogonal, but not necessarily a rotation matrix
		Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
		I(2, 2) = (U * V.transpose()).determinant();

		// Add the rotation matrix R to the list
		Eigen::Matrix3d R = U * I * V.transpose();
		all_RotationMatrices[i] = R;
	}
	return all_RotationMatrices;
}

/**
  * Compute b
  *
  * Computes the right hand side of the ecuation or the b vector using the formula from the paper
  *
  * @param rotationMatrices: vector with the rotation matrices of the vertices
  * @return MatrixXd b
  */
std::shared_ptr <Eigen::MatrixXd> ArapAlgorithm::compute_b(std::vector<Eigen::Matrix3d> rotationMatrices) {
	// Initialize the matrix b
	std::shared_ptr <Eigen::MatrixXd> b = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(verticesPreDeformation->rows(), 3));

	auto undeformedVertices = verticesPreDeformation;
	auto& neighborhood = neighbors;
	auto& weightMatrix = weight_Matrix;
	auto& fixedVertices = fixed_Vertices;

#ifdef OMP
#pragma omp parallel for default(none) \
            shared(undeformedVertices, neighborhood, weightMatrix, fixedVertices, rotationMatrices, b)
#endif

	// iterate over the vertices
	for (int i = 0; i < undeformedVertices->rows(); i++) {
		// create the row for the matrix b
		Eigen::Vector3d b_row = Eigen::Vector3d(0.0, 0.0, 0.0);
		// Don't forget to take into account fixed vertices!!!
		if (i == moving_Vertex) { // check to see if it is the vertex that is moving (saved in arap)
			b_row = movingVertex_Position; // set the value as the moving position also saved in arap
		}
		// try to see if the vertex is a fixed vertex by looking in the fixed Vertices matrix
		else if(std::find(fixedVertices.begin(), fixedVertices.end(), i) != fixedVertices.end()){
			b_row = undeformedVertices->row(i);
			// if it is fixed it needs to stay the same so take the pre deformed vertex position
		} // the vertex is not moving and is not fixed
		else { // vertex needs to be deformed
			for (int j : neighborhood[i]) { // Iterate over the neighbors
				// follow formula in paper: 1/2 * wij * (Ri + Rj)*(pi - pj)
				// compiler won't let me put it that way 
				b_row += 0.5 * weightMatrix(i,j) * (undeformedVertices->row(i) - undeformedVertices->row(j)) * 
					(rotationMatrices[i] + rotationMatrices[j]);
			}
		}
		b->row(i) = b_row;
	}
	return b;
}

/**
  * Calculate the rigidity energy
  *
  * Calculate rigidity energy of each iteration using the formula from paper
  * E(ci, C'i) = sum over neighbors ( wij * ||(p'i - p'j) - Ri*(pi - pj) ||^2)
  *
  * @param changedVertices: deformed vertices 
  * @param rotationMatrices: rotation matrices of the vertices
  * @return rigidity energy
  */
double ArapAlgorithm::calculateEnergy(std::shared_ptr<Eigen::MatrixXd> changedVertices, std::vector<Eigen::Matrix3d> rotationMatrices) {
	double rig_energy = 0.0; // rigidity energy
	auto undeformedVertices = verticesPreDeformation;
	auto& neighborhood = neighbors;
	auto& weightMatrix = weight_Matrix;

	double energyCell = 0.0; // energy per cell

#ifdef OMP
#pragma omp parallel for default(none) \
            shared(undeformedVertices, neighborhood, weightMatrix, changedVertices, rotationMatrices) \
            reduction(+: rig_energy) // makes it possible to do addition in parallel for the energy
#endif

	// iterate over vertices pre deformation (pi)
	for (int i = 0; i < undeformedVertices->rows(); i++) {
		energyCell = 0.0; // reset each time we switch vertex
		// E(ci, C'i) = sum over neighbors ( wij * ||(p'i - p'j) - Ri*(pi - pj) ||^2)
		for (int j : neighborhood[i]) { // Iterate over the neighbors
			// (p'i - p'j)
			Eigen::Vector3d deformed_edge = changedVertices->row(i) - changedVertices->row(j);
			// (pi - pj)
			Eigen::Vector3d undeformed_edge = undeformedVertices->row(i) - undeformedVertices->row(j);
			// solve equation for each cell
			energyCell += weightMatrix(i, j) * (deformed_edge - rotationMatrices[i] * undeformed_edge).squaredNorm();
		}
		// after finding the energy for each cell, add to the total energy
		rig_energy += weightMatrix(i, i) * energyCell; // reset each time we switch vertex
	}

	return rig_energy;
}

/**
  * Compute the arap algorithm
  *
  * Executes the arap algorithm following the main paper
  *
  * @param vertices: vertices in the mesh
  * @return MatrixXd new vertices
  */
std::shared_ptr<Eigen::MatrixXd> ArapAlgorithm::doArap(std::shared_ptr<Eigen::MatrixXd> vertices) {
	// save the previous 

	// initial guess
	solver.compute(system_Matrix_L.sparseView());
	std::shared_ptr<Eigen::MatrixXd> new_vertices = std::make_shared<Eigen::MatrixXd>(solver.solve(system_Matrix_L * *vertices));

	auto prev_rig_energy = DBL_MAX; 
	double actual_rig_energy = 0;

	// To test the times
	 // Start the timer
	const std::chrono::time_point<std::chrono::system_clock> t0 = std::chrono::system_clock::now();

	// ITERATE
	int k = 0;
	do {
		// calculate rotation matrices
		std::vector<Eigen::Matrix3d> allRotationMatrices = calculateRotationMatrix(new_vertices); // this has to be the deformed vertices

		// Compute b once we have the rotation matrices
		std::shared_ptr < Eigen::MatrixXd> b = compute_b(allRotationMatrices);

		// update mesh
		verticesPreDeformation = new_vertices;

		//find optimal mesh
		solver.compute(system_Matrix_L.sparseView());
		new_vertices = std::make_shared<Eigen::MatrixXd>(solver.solve(*b));

		// update energies
		actual_rig_energy = calculateEnergy(new_vertices, allRotationMatrices);

		k++;
		if (k >= 4 && abs(prev_rig_energy - actual_rig_energy) < 8e-1) 
		{
			DEBUG_STDOUT("Iteration %d: Energy threshold %.4f reached! Stopping early..\n", k, 8e-1);
			k = 10;
		}
		else {
			prev_rig_energy = actual_rig_energy;
		}
	} while (k < 10);

	// print the highest numbers
	DEBUG_STDOUT("Highest iteration %d: , Highest rigidity energy = %.4f\n", k, actual_rig_energy);

	// End the timer and print the duration
	const std::chrono::time_point<std::chrono::system_clock> t1 = std::chrono::system_clock::now();
	printf("\n\nTook %f seconds to deform..\n\n", std::chrono::duration<double>(t1 - t0).count());
	printf("Number of vertices: %d", new_vertices->rows());

	return new_vertices;
}