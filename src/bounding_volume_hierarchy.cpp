#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>

/* rootNodeHelper - Helper function to create the root node.
*  Inputs:
*   - pScene: pointer to the scene
*   - indexes: vector of indexes that is stored by the root node, containing the mesh and index within the mesh for each triangle in the scene
*   - x_low, x_high, y_low, y_high, z_low, z_high: variables that define the outer boundaries of the AABB.
*  Outputs: None. The input variables will be updated to fit the scene.
*/
void rootNodeHelper(Scene* pScene, std::vector<int> & indexes, float & x_low, float & x_high, float & y_low, float & y_high, float & z_low, float & z_high) {
    for (int i = 0; i < pScene->meshes.size(); ++i) {
        Mesh mesh = pScene->meshes[i];
        for (int j = 0; j < mesh.triangles.size(); ++j) {
            indexes.push_back(i); // Mesh in the scene
            indexes.push_back(j); // Triangle in mesh

            glm::uvec3 triangle = mesh.triangles[j];
            std::vector<Vertex> vertices { mesh.vertices[triangle.x], mesh.vertices[triangle.y], mesh.vertices[triangle.z] };
            for (int k = 0; k < 3; ++k) {
                x_low = fmin(x_low, vertices[k].position.x);
                x_high = fmax(x_high, vertices[k].position.x);
                y_low = fmin(y_low, vertices[k].position.y);
                y_high = fmax(y_high, vertices[k].position.y);
                z_low = fmin(z_low, vertices[k].position.z);
                z_high = fmax(z_high, vertices[k].position.z);
            }
        }
    }
}

/* mergeSortBy - perform merge sort on a vector, changing the order of the other in the process.
*  Inputs:
*   - by: vector of size n to be sorted
*   - A: vector of size n. If value x on position i in by gets assigned the new position j, the value at A[i] will also move to A[j]
*/
void mergeSortBy(std::vector<float> &by, std::vector<glm::uvec3> &A) {
    assert(by.size() == A.size());
}


void BoundingVolumeHierarchy::growBVH(int nodeIndex, int recursionDepth) {
    if (this->nodes[nodeIndex].indexes.size() == 2 || recursionDepth >= this->maxLevels) {
        /* Exit condition 1: the node contains vector entries for 1 mesh and 1 triangle (size of 2 entries)
         * Exit condition 2: the tree is at its maximum recursion depth.
         * In either case, the node is a leaf node.
         */ 
        if (recursionDepth > this->m_numLevels)
            this->m_numLevels = recursionDepth;

        this->m_numLeaves = this->m_numLeaves + 1;
        this->nodes[nodeIndex].isParent = false;
    } else {
        // For nodes in the middle of the tree ...
        Node thisNode = this->nodes[nodeIndex];

        int dimension = recursionDepth % 3;             // 0 = x, 1 = y, 2 = z
        std::vector<glm::uvec3> triangles;
        std::vector<std::vector<Vertex>> vertices;
        std::vector<float> means;

        // Extract the center positions of the triangles along the dimension of interest.
        for (int n = 0; n < this->nodes[nodeIndex].indexes.size() / 2; ++n) {
            Mesh thisMesh = this->m_pScene->meshes[thisNode.indexes[2 * n]];
            triangles[n] = thisMesh.triangles[thisNode.indexes[2*n+1]];
            vertices[n][0] = thisMesh.vertices[triangles[n].x];
            vertices[n][1] = thisMesh.vertices[triangles[n].y];
            vertices[n][2] = thisMesh.vertices[triangles[n].z];
            means[n] = (vertices[n][0].position[dimension] + vertices[n][1].position[dimension] + vertices[n][2].position[dimension]) / 3;
                        // !!! For debugging: the line above uses position[dimension] to access .x,.y,.z Is this supposed to work?
        }

        // Sort the triangles vector using means.
        mergeSortBy(means,triangles); // To be continued...
    }
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    // Create root node containing all triangles in the scene and an AABB for the entire scene, set to false to indicate it's a leaf node.
    std::vector<int> indexes;
    float x_low = std::numeric_limits<float>::infinity(), y_low=x_low, z_low=x_low;
    float x_high = - std::numeric_limits<float>::infinity(), y_high=x_high, z_high=x_high;
    rootNodeHelper(m_pScene, indexes, x_low, x_high, y_low, y_high, z_low, z_high);
    Node root(true, x_low, x_high, y_low, y_high, z_low, z_high, indexes);

    // Initialize the 0th layer of the BVH.
    nodes = { root };
    m_numLevels = 0;
    m_numLeaves = 1;

    // Recursively create the BVH.
    growBVH(0,0);
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return 1;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return 1;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
}


// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}


// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        bool hit = false;
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;
                    hit = true;
                }
            }
        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        return false;
    }
}