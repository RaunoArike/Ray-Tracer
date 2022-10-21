#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>

#include <iostream>

/* findBounds - Helper function to find boundaries of an AABB in the scene.
*  Inputs:
*   - pScene - pointer to the scene 
*   - indexes - indexes array according to format [mesh0, triangle0, mesh1, triangle1, ...] as defined in the Node struct.
*  Outputs:
*   - bounds - 2D vector that holds the boundaries of the AABB.
*/
std::vector<std::vector<float>> findBounds(Scene* pScene, std::vector<int> indexes) {
    assert(!(indexes.size() == 0));

    // Definition, setting bounds to +- infinity for default.
    std::vector<std::vector<float>> bounds = { { std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity() },
        { std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity() },
        { std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity() } };

    // Adjust the boundaries by iterating over the triangles.
    int i = 0;
    int meshNo = *indexes.begin();                              // Current mesh number and copy of it.
    Mesh mesh = pScene->meshes[meshNo];
    for (auto it = indexes.begin(); it != indexes.end(); ++it,++i) {
        // Reload the mesh data type if it's a new mesh. Prevents copying a lot of data in larger meshes.
        if (*it != meshNo) {
            meshNo = *it;
            mesh = pScene->meshes[meshNo];
        }

        // Find the triangle within the mesh.
        ++it;
        glm::uvec3 triangle = mesh.triangles[*it];

        // Find the vertices that compose this triangle. Adjust the bounds according to the individual x,y,z values of this AABB.
        for (int k = 0; k < 3; k++) {
            Vertex v = mesh.vertices[triangle[k]];      

            float x = v.position.x;
            float y = v.position.y;
            float z = v.position.z;

            bounds[0][0] = fmin(bounds[0][0], x);
            bounds[0][1] = fmax(bounds[0][1], x);
            bounds[1][0] = fmin(bounds[1][0], y);
            bounds[1][1] = fmax(bounds[1][1], y);
            bounds[2][0] = fmin(bounds[2][0], z);
            bounds[2][1] = fmax(bounds[2][1], z);
        }
    }
    return bounds;

    /* Remaining problem with this function: it iterates over the indexes (=2*n_triangles) instead of the vertices, even though the 
    *  vertices decide the final x,y,z values of the AABB, wasting time because the amount of triangles is larger than the amount of
    *  vertices. 
    */
}

/* rootNodeHelper - Helper function to create the root node.
*  Inputs:
*   - pScene: pointer to the scene
*   - indexes: vector of indexes that is stored by the root node, containing the mesh and index within the mesh for each triangle in the scene
*   - bounds: 2D vector that holds the boundaries of the AABB.
*  Outputs: None. The input variables will be updated to fit the scene.
*/
void rootNodeHelper(Scene* pScene, std::vector<int> & indexes, std::vector<std::vector<float>> & bounds) {
    for (int i = 0; i < pScene->meshes.size(); ++i) {
        Mesh mesh = pScene->meshes[i];
        for (int j = 0; j < mesh.triangles.size(); ++j) {
            indexes.push_back(i); // Mesh in the scene
            indexes.push_back(j); // Triangle in mesh
        }
    }
    bounds = findBounds(pScene, indexes);
}

/* merge - performs merging operation for merge sort.
*  Inputs:
*   - by - vector size n that indicates the order
*   - byL, byR - vectors that are to be merged into by
*   - A - vector size 2n that is sorted according to the order in by
*   - aL, aR - vectors that are merged into A
*  Outputs:
*   None. Variables are passed by reference.
*/
void merge(std::vector<float> &by, std::vector<float>& byL, std::vector<float>& byR, std::vector<int> &A, std::vector<int>& aL, std::vector<int>& aR) {
    // Scan through both vectors after modification, selecting the smallest in order every time, changing the order of A accordingly. 
    int i = 0, j = 0;
    while (i < byL.size() && j < byR.size()) {
        if (byR[j] < byL[i]) {
            by.push_back(byR[j]);
            A.push_back(aR[2 * j]);
            A.push_back(aR[2 * j + 1]);
            j++;
        } else { 
            by.push_back(byL[i]);
            A.push_back(aL[2 * i]);
            A.push_back(aL[2 * i + 1]);
            i++;
        }
    }
    while (i >= byL.size() && j < byR.size()) {
        by.push_back(byR[j]);
        A.push_back(aR[2 * j]);
        A.push_back(aR[2 * j+1]);
        j++;
    }
    while (i < byL.size() && j >= byR.size()) {
        by.push_back(byL[i]);
        A.push_back(aL[2 * i]);
        A.push_back(aL[2 * i+1]);
        i++;
    }
}

/* mergeSortBy - perform merge sort in ascending direction on distance vector, changing the order of the indexes vector in the process.
*  Inputs:
*   - by: vector of size n to be sorted
*   - A: vector of size 2n. If value x on position i in 'by' gets assigned the new position j, the values at A[2i,2i+1] will also move to A[2j,2j+1]
*  Outputs:
*   - None. The input vectors are passed by reference and are modified accordingly by the function.
*/
void mergeSortBy(std::vector<float> &by, std::vector<int> &A) {
    assert(by.size() == A.size()/2);

    if (!(by.size() <= 1)) {                // The function is not executed if array size <= 1 (=exit condition)
        std::vector<float> byL;             // left split &by
        std::vector<float> byR;             // right split &by
        std::vector<int> aL;                // left split &A
        std::vector<int> aR;                // right split &A

        // Put the means and indexes of the vectors into the left and right subvectors.
        // Erase by and A in the process
        int half = by.size() / 2;
        for (int i = 0; i < by.size(); ++i)
        {
            if (i < half) {
                byL.push_back(by[i]);
                aL.push_back(A[2 * i]);
                aL.push_back(A[2 * i + 1]);
            } else {
                byR.push_back(by[i]);
                aR.push_back(A[2 * i]);
                aR.push_back(A[2 * i + 1]);
            }
        }
        // Empty by and A 
        by.clear();
        A.clear();

        // Recursive call
        mergeSortBy(byL, aL);
        mergeSortBy(byR, aR);

        // Merge
        merge(by, byL, byR, A, aL, aR);
    } 
}

/* growBVH - uses recursion to grow the bounding volume hierarchy
*  Inputs: 
*   - nodeIndex - index of the node in the nodes vector in the BVH
*   - recursionDepth - depth of the current recursion step
*  Outputs:
*   - None. Modifies the attributes of the BVH to shape it 
*/
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
        Node thisNode = this->nodes[nodeIndex];          // Pass by reference !
        std::vector<int> indexes = thisNode.indexes;

        int dimension = recursionDepth % 3;             // 0 = x, 1 = y, 2 = z
        std::vector<glm::uvec3> triangles;
        std::vector<std::vector<Vertex>> vertices;
        std::vector<float> means;

        // Extract the center positions of the triangles along the dimension of interest.
        int meshNo = indexes[0];
        Mesh thisMesh = this->m_pScene->meshes[meshNo];
        for (int n = 0; n < indexes.size() / 2; ++n) {
            if (indexes[2*n] != meshNo)                                         // Load new mesh iff the index at 2n points to a new one. Spares a lot of time copying data.
                Mesh thisMesh = this->m_pScene->meshes[indexes[2 * n]];

            triangles.push_back(thisMesh.triangles[indexes[2*n+1]]);
            vertices.push_back({ thisMesh.vertices[triangles[n].x], 
                                 thisMesh.vertices[triangles[n].y], 
                                 thisMesh.vertices[triangles[n].z] });
            means.push_back((vertices[n][0].position[dimension] + vertices[n][1].position[dimension] + vertices[n][2].position[dimension]) / 3);
                        // !!! For debugging: the line above uses position[dimension] to access .x,.y,.z Is this supposed to work?
        }

        // Sort the triangles vector using means.
        mergeSortBy(means,indexes);

        // Find the median in the middle of the sorted list.
        int medianIndex = means.size() / 2 - 1;
        
        // Split up the indexes vector for the left and right child node
        std::vector<int> indexesL, indexesR;
        int i = 0;
        while (i <= medianIndex) {                      // Left node holds all lower triangles up to and including the median.
            indexesL.push_back(indexes[2 * i]);
            indexesL.push_back(indexes[2 * i + 1]);
            ++i;
        }
        while (i < means.size()) {                      // Right node takes the triangles higher than the median.
            indexesR.push_back(indexes[2 * i]);
            indexesR.push_back(indexes[2 * i + 1]);
            ++i;
        }

        // Set the other variables for the child nodes.
        std::vector<std::vector<float>> boundsL = findBounds(this->m_pScene, indexesL);
        std::vector<std::vector<float>> boundsR = findBounds(this->m_pScene, indexesR);

        Node left = { true, boundsL, indexesL };
        Node right = { true, boundsR, indexesR };

        // Put the nodes in the node array of the BVH.
        int leftIndex = this->nodes.size();
        this->nodes[nodeIndex].indexes = { leftIndex, leftIndex + 1 };           // thisNode was passed by reference
        this->nodes.push_back(left);
        this->nodes.push_back(right);

        // Recursive call for both child nodes.
        this->growBVH(leftIndex, recursionDepth + 1);       // Left node
        this->growBVH(leftIndex+1, recursionDepth + 1);     // Right node
    }
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    // Create root node containing all triangles in the scene and an AABB for the entire scene, set to false to indicate it's a leaf node.
    std::vector<int> indexes;
    std::vector<std::vector<float>> bounds;
    rootNodeHelper(m_pScene, indexes, bounds);
    Node root(true, bounds, indexes);

    // Initialize the 0th layer of the BVH.
    this->nodes = { root };
    this->m_numLevels = 0;
    this->m_numLeaves = 1;

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