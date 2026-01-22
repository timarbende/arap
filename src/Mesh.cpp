/**
 * Project: Interactive ARAP
 * File:    Mesh.cpp
 * Authors: Ankur Deria, Michael Dey, Bendeguz Timar, Andrea Solanas de Vicente
 */

#include "Mesh.h"

int jump_lerp_steps = 50;
//float distance = .2;
float distance = .05;
static Eigen::Vector3f jumpVector = Eigen::Vector3f(distance, distance, distance) / jump_lerp_steps;
int jump_lerp=-1;

/**
  * Mesh constructor 3
  *
  * Constructs the mesh by an off file with the vertices and faces
  * 
  * @param path: file name of the location where the mesh is stored
  * @return Mesh instance
  */
Mesh::Mesh(const std::string& path) {
	// Load a mesh in OFF format
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	if (!igl::readOFF(path, V, F))
	{
		DEBUG_STDERR("Failed to read file: ");
		DEBUG_STDERR(path);
	}
	vertices = std::make_shared<Eigen::MatrixXd>(std::move(V));
	faces = std::make_shared<Eigen::MatrixXi>(std::move(F));
	// Initialize white colors for all the faces
	colors = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Constant(faces->rows(), 3, 1));
}

/**
  * Find closest vertex to the selected face
  *
  * Calculates through maximum coefficient which is the vertex closest to the selected face
  *
  * @param face_id: id of the face that has been selected
  * @param position: the barycentric position of the mouse
  * @return id of the vertex closest to the clicked position
  */
int Mesh::findClosestVertexToSelectedFace(int face_id, const Eigen::Vector3f& position) {
	// the vertex with the highest coefficient is the closest 
	// ( alpha * vertex1 + beta * vertex2 + gamma * vertex3 = position of face)
	int max_n;
	position.maxCoeff(&max_n);
	int sol = faces->row(face_id)(max_n);
	return sol;
}

/**
  * Convert camera position to world position
  *
  * Converts the position in camera space to world space 
  *
  * @param viewer: libigl viewer
  * @param vertex_id: the vertex id of which we have to convert the position
  * @return vector of position
  */
Eigen::Vector3f Mesh::convertCameraToWorldPosition(igl::opengl::glfw::Viewer& viewer, int vertex_id) {
	// otain the mouse position
	Eigen::Vector2f mousePosition = GetMousePosition(viewer);
	// obtain the vertex position using the vertex id
	Eigen::Vector3f vertexPosition = { (float)vertices->row(vertex_id).x(), (float)vertices->row(vertex_id).y(), (float)vertices->row(vertex_id).z() };
	// calculate the projection using igl: from 3d to 2d
	Eigen::Vector3f projection = igl::project(vertexPosition, viewer.core().view, viewer.core().proj, viewer.core().viewport);
	// obtain the world positon using the mouse position
	Eigen::Vector3f worldPosition = igl::unproject(Eigen::Vector3f(mousePosition.x(), mousePosition.y(), projection.z()),
		viewer.core().view, viewer.core().proj, viewer.core().viewport);

	return worldPosition;
}

/**
  * Use Arap algorithm
  *
  * Runs the ARAP algorithm using the input from the user
  *
  * @param viewer: libigl viewer
  * @return void
  */
void Mesh::useARAP(igl::opengl::glfw::Viewer& viewer){
	// collect the fixed faces for arap
	fixed_faces.clear();
	fixed_faces.reserve(fixed_faces_selected.size());
	for (auto entry : fixed_faces_selected) {
		if (entry.second) { // If face is selected
			fixed_faces.push_back(entry.first);
		}
	}
	
	// Using class ArapAlgorithm
	// saves in that class the new vertex and its position, then obtain the fixed vertices of the faces, update the system matrix L
	aa.updateMovingVertex(moving_vertex, convertCameraToWorldPosition(viewer, moving_vertex), faces, fixed_faces);
	if (jump_lerp != -1) {
		Eigen::Vector3f vertexPosition = { (float)vertices->row(moving_vertex).x(), (float)vertices->row(moving_vertex).y(), (float)vertices->row(moving_vertex).z() };
		aa.updateMovingVertex(moving_vertex, vertexPosition + jumpVector, faces, fixed_faces);
	}

	// Compute deformation
	std::shared_ptr<Eigen::MatrixXd> deformedVertices = std::move(aa.doArap(vertices));
	vertices = std::make_shared< Eigen::MatrixXd>(*deformedVertices);

	// Change the mesh in the viewer
	viewer.data().compute_normals();
	viewer.data().set_mesh(*vertices, *faces);
}

/**
  * Handle selection in mesh
  *
  * Handles the selection of a face in the displayed mesh
  *
  * @param viewer: libigl viewer
  * @return true if the face was registered correctly, false if it wasn't
  */
bool Mesh::handleSelection(igl::opengl::glfw::Viewer& viewer){
	// reset the face id and position
	faceId = -1;
	barycentricPosition = Eigen::Vector3f(0, 0, 0);

	if (isArapOn){ // if the algorithm is working
		if (moving_vertex == -1){   // If no moving vertex has been selected yet
			if (igl::unproject_onto_mesh(GetMousePosition(viewer), viewer.core().view, viewer.core().proj, viewer.core().viewport, *vertices, *faces, faceId, barycentricPosition))	{
				// using the mouse psoition and the face id, find the closest vertex (this is what arap algorithm uses)
				moving_vertex = findClosestVertexToSelectedFace(faceId, barycentricPosition);

				// Pre compute for ARAP algorithm
				aa.preComputeDeformation(vertices, faces);
				return true;
			}
		}
	}
	else { // arap is not running
		if (igl::unproject_onto_mesh(GetMousePosition(viewer), viewer.core().view, viewer.core().proj, viewer.core().viewport,
			*vertices, *faces, faceId, barycentricPosition)) {

			// Store the selections of the faces
			if (fixed_faces_selected.contains(faceId)) { // if it was already selected change it
				fixed_faces_selected[faceId] = !fixed_faces_selected[faceId];
			}
			else { // if not selected put it as fixed
				fixed_faces_selected[faceId] = true; 
			}
			
			// Set the color for the selected face -> red
			colors->row(faceId) << 1, !fixed_faces_selected[faceId], !fixed_faces_selected[faceId];

			// Paint
			this->viewer.data().set_colors(*colors);
			return true;
		}
	}
	return false;
}

/**
  * Get the position of the mouse
  *
  * Obtain the position of the mouse in the viewer
  *
  * @param viewer: libigl viewer
  * @return vector of position (2D)
  */
Eigen::Vector2f Mesh::GetMousePosition(igl::opengl::glfw::Viewer& viewer) {
	return Eigen::Vector2f{ viewer.current_mouse_x, viewer.core().viewport(3) - (float)viewer.current_mouse_y };
}

/**
  * Handle mouse clicked 
  *
  * Handles the mouse being clicked or mouse down
  *
  * @param viewer: libigl viewer
  * @param button: button that was pressed
  * @param modifier:
  * @return true if the event was handled correctly, false if it wasn't
  */
bool Mesh::handleMouseDownEvent(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
	hasMeshBeenSelected = handleSelection(viewer);

	return hasMeshBeenSelected;
}

/**
  * Handle mouse being unclicked
  *
  * Handles the mouse not being clicked anymore or mouse up
  *
  * @return void
  */
void Mesh::handleMouseReleaseEvent() 
{
	hasMeshBeenSelected = false; // set to false that the mesh has been selected
	moving_vertex = -1; // Reset the selection of moving vertex
}

/**
  * Handle mouse moving
  *
  * Handles the mouse moving once it has clicked on a face
  *
  * @param viewer: libigl viewer
  * @return true if the event was handled correctly, false if it wasn't
  */
bool Mesh::handleMouseMoveEvent(igl::opengl::glfw::Viewer& viewer) {
	if (hasMeshBeenSelected) {
		handleSelection(viewer); // Selecting vertices for deformation
		if (isArapOn) {
			useARAP(viewer); // Do ARAP, Do the deformation 
		}
		return true;
	}
	return false;
}

/**
  * Handle key being pressed
  *
  * Handles a key from the keyboard being presed down (important keys in this are A and B)
  *
  * @param keyPressed: character of the key that was pressed
  * @return true if the event was handled correctly, false if it wasn't
  */
bool Mesh::handleKeyDownEvent(unsigned char keyPressed) {
	if (keyPressed == 'A') {
		isArapOn = !isArapOn; // to turn on Arap if it wasn't on
		for (auto& a : fixed_faces_selected) {
			if (a.second) { // if the face was selected
				colors->row(a.first) << !isArapOn, isArapOn, 0; // change color of faces -> green
			}
		}
		// set the colors now
		viewer.data().set_colors(*colors);
		return true;
	} // if they click another key
	else if (keyPressed == 'B') { // reset button
		isArapOn = false; // no more arap
		// remove the fixed faces
		fixed_faces_selected.clear();
		// put the colors in white again
		colors = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Constant(faces->rows(), 3, 1));
		viewer.data().set_colors(*colors);
		return true;
	}
	else if (keyPressed == 'C') { // jumpMovingVector
		if (moving_vertex != -1) {
			for (jump_lerp = jump_lerp; jump_lerp < jump_lerp_steps; jump_lerp++) {
				useARAP(viewer);
			}
			jump_lerp = -1;
		}
		return true;
	}
	// in the case of other letters
	isArapOn = false;
	return false;
}

/**
  * Launch the libigl viewer
  *
  * Launches the viewer for the mesh
  *
  * @return void
  */
void Mesh::launchViewer(){
	// handle the mouse being clicked
	viewer.callback_mouse_down = [this](igl::opengl::glfw::Viewer& viewer, int button, int modifier)->bool {return Mesh::handleMouseDownEvent(viewer, button, modifier); };
	// handle the mouse moving while being clicked
	viewer.callback_mouse_move = [this](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {return Mesh::handleMouseMoveEvent(viewer); };
	// handle a key being pressed
	viewer.callback_key_down = [this](igl::opengl::glfw::Viewer& viewer, unsigned char keyPressed, int) -> bool { return Mesh::handleKeyDownEvent(keyPressed); };
	// handle the mouse not being clicked anymore
	viewer.callback_mouse_up = [this](igl::opengl::glfw::Viewer& viewer, int button, int modifier)->bool {
		Mesh::handleMouseReleaseEvent();
		return true;
	};

	// Plot the mesh
	viewer.data().set_mesh(*vertices, *faces);
	viewer.data().set_colors(*colors);
	// Launch the glfw viewer
	viewer.launch();
}


