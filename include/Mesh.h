/**
 * Project: Interactive ARAP
 * File:    Mesh.h
 * Authors: Ankur Deria, Michael Dey, Bendeguz Timar, Andrea Solanas de Vicente
 */

#pragma once

#include <cstdlib>;
#include <memory>;
#include <Eigen/Dense>;
#include <ArapAlgorithm.h>;
#include <igl/readOFF.h>;
#include <igl/opengl/glfw/Viewer.h>;
#include <igl/unproject_onto_mesh.h>;
#include <igl/project.h>;
#include <igl/unproject.h>;

 /*
  * Class that implements a Mesh
  */
class Mesh {
public:
	// Constructor 1
	Mesh() = default;
	// Constructor 2
	Mesh(Eigen::MatrixXd&& V, Eigen::MatrixXi&& F) : vertices(std::make_shared<Eigen::MatrixXd>(std::move(V))), faces(std::make_shared<Eigen::MatrixXi>(std::move(F))) {}
	// Constructor 3
	Mesh(const std::string& path);

	// launch the viewer
	void launchViewer();
	
private:

	// VARIABLES 
	std::shared_ptr<Eigen::MatrixXd> vertices; // Vertices from the mesh
	std::shared_ptr<Eigen::MatrixXi> faces; // Faces from the mesh
	std::shared_ptr<Eigen::MatrixXd> colors; // Colors from the mesh
	igl::opengl::glfw::Viewer viewer{}; // Create a viewer object and set the mesh to be displayed
	std::vector<int> fixed_faces; // ids of the fixed faces in the mesh, then will become vertices
	std::map<int, bool> fixed_faces_selected; // Selected faces (fixed) -> id of face and true if they are selected
	int faceId; // face ide selected, used when user selects all the faces
	ArapAlgorithm aa = ArapAlgorithm(); // Arap instance for the algorithm
	int moving_vertex = -1; // Selected moving vertex to be used to perform ARAP, set to -1 when nothing is selected
	Eigen::Vector3f barycentricPosition; // position of the moving vertex (braycentric)
	bool hasMeshBeenSelected = false; // If clicked on a point on the mesh it is true, else false
	bool isArapOn = false; // If ARAP is running it's true, else false
	bool isMouseDown; // if mouse is pressed down

	// FUNCTIONS
	// Function to execute the Arap algorithm once the faces are selected and the moing vertex has been selected as well
	void useARAP(igl::opengl::glfw::Viewer& viewer);
	// Returns the mouse position
	static Eigen::Vector2f GetMousePosition(igl::opengl::glfw::Viewer& viewer);
	// Handles the selection (finds and stores the selected face and the closest vertex to the selection)
	bool handleSelection(igl::opengl::glfw::Viewer& viewer);
	// Event listeners
	bool handleKeyDownEvent(unsigned char keyPressed);
	void handleMouseReleaseEvent();
	bool handleMouseMoveEvent(igl::opengl::glfw::Viewer& viewer);
	bool handleMouseDownEvent(igl::opengl::glfw::Viewer& viewer, int button, int modifier);
	// Finds the closest vertex index to a selected face
	int findClosestVertexToSelectedFace(int face_id, const Eigen::Vector3f& position);
	// Converts the camera position of the vertex to a world position
	Eigen::Vector3f convertCameraToWorldPosition(igl::opengl::glfw::Viewer& viewer, int vertex_id);
};

