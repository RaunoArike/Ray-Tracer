#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>
#include <screen.cpp>
#include <queue>


/*** Helper functions ***/


/* findBounds - Helper function to find boundaries of an AABB in the scene.
*  Inputs:
*   - pScene - pointer to the scene 
*   - indexes - indexes array according to format [mesh0, triangle0, mesh1, triangle1, ...] as defined in the Node struct.
*  Outputs:
*   - bounds - 2D vector that holds the boundaries of the AABB.
*/
std::vector<std::vector<float>> findBounds(Scene* pScene, std::vector<int> &indexes) {
    // Definition, setting bounds to +- infinity for default.
    std::vector<std::vector<float>> bounds = { { std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity() },
        { std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity() },
        { std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity() } };

    // Adjust the boundaries by iterating over the triangles.
    int i = 0;
    for (auto it = indexes.begin(); it != indexes.end(); ++it) {
        int meshNo = *it;          // Register index of the mesh in the scene before moving to triangle #
        ++it;

        glm::uvec3 triangle = pScene->meshes[meshNo].triangles[*it];

        // Find the vertices that compose this triangle. Adjust the bounds according to the individual x,y,z values of this AABB.
        for (int k = 0; k < 3; k++) {
            Vertex v = pScene->meshes[meshNo].vertices[triangle[k]];

            bounds[0][0] = fmin(bounds[0][0], v.position.x);
            bounds[0][1] = fmax(bounds[0][1], v.position.x);
            bounds[1][0] = fmin(bounds[1][0], v.position.y);
            bounds[1][1] = fmax(bounds[1][1], v.position.y);
            bounds[2][0] = fmin(bounds[2][0], v.position.z);
            bounds[2][1] = fmax(bounds[2][1], v.position.z);
        }
    }
    return bounds;

    /* Remaining problem with this function: it iterates over the indexes (=2*n_triangles) instead of the vertices, even though the 
    *  vertices decide the final x,y,z values of the AABB, wasting time because the amount of triangles is larger than the amount of
    *  vertices. (and we check way more vertices than needed, because triangles may point to the same vertices.)
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
        for (int j = 0; j < pScene->meshes[i].triangles.size(); ++j) {
            indexes.push_back(i); // Mesh in the scene
            indexes.push_back(j); // Triangle in mesh
        }
    }
    bounds = findBounds(pScene, indexes);
}

/* findMeans - create vector means that holds the means of triangles in pScene at each index along a specific dimension.
 *  Inputs:
 *   - pScene - pointer to the scene
 *   - indexes -
 */
std::vector<float> findMeans(Scene* pScene, std::vector<int>& indexes, int dimension)
{
    std::vector<float> means;
    for (int n = 0; n < indexes.size() / 2; ++n) {
        glm::uvec3 triangle = pScene->meshes[indexes[2 * n]].triangles[indexes[2 * n + 1]];
        std::vector<Vertex> vertices = { pScene->meshes[indexes[2 * n]].vertices[triangle.x],
            pScene->meshes[indexes[2 * n]].vertices[triangle.y],
            pScene->meshes[indexes[2 * n]].vertices[triangle.z] };
        means.push_back((vertices[0].position[dimension] + vertices[1].position[dimension] + vertices[2].position[dimension]) / 3);
    }
    return means;
}

/* merge - performs merging operation for merge sort.
 *  Inputs:
 *   - by - vector size n that indicates the order
 *   - A - vector size 2n that is sorted according to the order in by
 *   - left - index of the leftmost variable under consideration in by
 *   - mid - index of the middle variable under consideration in by
 *   - right - index of the rightmost variable under consideration in by
 *  Outputs:
 *   None. Variables are passed by reference and modified accordingly. When a variable at by[i] is moved to by[j], then variables A[2i] and A[2i+1] are moved to A[2j] and A[2j+1] respectively.
 */
void mergeBy(std::vector<float>& by, std::vector<int> &A, int const left, int const mid, int const right)
{
    std::vector<float> byTmp(right-left+1);       // Preallocate space for temporary vectors to save time.
    std::vector<int> aTmp(2 * (right-left+1));
    int i = left, j = mid + 1, k = 0;

    while (i <= mid && j <= right) {
        if (by[i] <= by[j]) {
            byTmp[k] = by[i];
            aTmp[2*k] = A[2*i];
            aTmp[2*k+1] = A[2 * i + 1];
            ++i;
            ++k;
        } else {
            byTmp[k] = by[j];
            aTmp[2*k] = A[2 * j];
            aTmp[2*k+1] = A[2 * j + 1];
            ++j;
            ++k;
        }
    }
    while (i <= mid) {
        byTmp[k] = by[i];
        aTmp[2*k] = A[2 * i];
        aTmp[2*k+1] = A[2 * i + 1];
        ++i;
        ++k;
    }
    while (j <= right) {
        byTmp[k] = by[j];
        aTmp[2*k] = A[2 * j];
        aTmp[2*k+1] = A[2 * j + 1];
        ++j;
        ++k;
    }
    for (int i = left; i <= right; ++i) {
        by[i] = byTmp[i - left];
        A[2 * i] = aTmp[2 * (i - left)];
        A[2 * i + 1] = aTmp[2 * (i - left) + 1];
    }

}

/* mergeSortBy - perform merge sort in ascending direction on distance vector, changing the order of the indexes vector in the process.
 *  Inputs:
 *   - by: vector of size n to be sorted
 *   - A: vector of size 2n. If value x on position i in 'by' gets assigned the new position j, the values at A[2i,2i+1] will also move to A[2j,2j+1]
 *  Outputs:
 *   - None. The input vectors are passed by reference and are modified accordingly by the function.
 * 
 *  Implementation of mergeSortBy() and mergeBy() are based on the implementation on: https://slaystudy.com/c-merge-sort-vector/ 
 */
void mergeSortBy(std::vector<float>& by, std::vector<int> &A, int const begin, int const end) {
    assert(by.size() == A.size() / 2);

    if (begin >= end) {
        return;
    } else {
        int mid = begin + (end - begin) / 2;
        mergeSortBy(by, A, begin, mid);
        mergeSortBy(by, A, mid + 1, end);
        mergeBy(by, A, begin, mid, end);
    }
}


/*** BoundingVolumeHierarchy class functions. ***/


/* growBVH - uses recursion to grow the bounding volume hierarchy
*  Inputs: 
*   - nodeIndex - index of the node in the nodes vector in the BVH
*   - recursionDepth - depth of the current recursion step
*  Outputs:
*   - None. Modifies the attributes of the BVH to shape it 
*/
void BoundingVolumeHierarchy::growBVH(int nodeIndex, int recursionDepth) {
    /* Exit condition 1: the node contains vector entries for 1 mesh and 1 triangle (size of 2 entries)
     * Exit condition 2: the tree is at its maximum recursion depth.
     * In either case, the node is a leaf node.
     */ 
    if (this->nodes[nodeIndex].indexes.size() == 2 || recursionDepth >= this->maxLevels) {
        if (recursionDepth > this->m_numLevels)
            this->m_numLevels = recursionDepth;

        this->m_numLeaves = this->m_numLeaves + 1;
        this->nodes[nodeIndex].isParent = false;
    } else {
        // Get indexes of this node.
        std::vector<int> indexes = this->nodes[nodeIndex].indexes;

        // Find means of the triangles among these indexes along the dimension of interest.
        std::vector<float> means = findMeans(this->m_pScene, indexes, recursionDepth % 3);

        // Sort the triangles vector by means.
        mergeSortBy(means,indexes,0,means.size()-1);

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

        Node left = { true, findBounds(this->m_pScene, indexesL), indexesL };
        Node right = { true, findBounds(this->m_pScene, indexesR), indexesR };

        // Put the nodes at the end of the node array of the BVH.
        int leftIndex = this->nodes.size();
        this->nodes[nodeIndex].indexes = { leftIndex, leftIndex + 1 };           // thisNode was passed by reference
        this->nodes.push_back(left);
        this->nodes.push_back(right);

        // Recursive call for both child nodes.
        this->growBVH(leftIndex, recursionDepth + 1);       // Left node
        this->growBVH(leftIndex+1, recursionDepth + 1);     // Right node
    }
}

/* BoundingVolumeHierarchy - class constructor
*  Outputs:
*   - Creates an instance of a BHV class.
*/
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
    this->m_numLeaves = 0;

    // Recursively create the BVH.
    growBVH(0,0);
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return this->m_numLevels;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return this->m_numLeaves;
}

void BoundingVolumeHierarchy::debugDrawLevel(int level, Node node) {
    if (level == 0 || !node.isParent) {
        AxisAlignedBox aabb { glm::vec3(node.bounds[0][0], node.bounds[1][0], node.bounds[2][0]),
            glm::vec3(node.bounds[0][1], node.bounds[1][1], node.bounds[2][1]) };

        // Draw the AABB as a transparent green box.
        // drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

        // Draw the AABB as a (white) wireframe box.
        // drawAABB(aabb, DrawMode::Wireframe);
        drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    } else {
        --level;
        debugDrawLevel(level, this->nodes[node.indexes[0]]);
        debugDrawLevel(level, this->nodes[node.indexes[1]]);
    }
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    assert(level <= this->maxLevels);
    debugDrawLevel(level, this->nodes[0]);
}


// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    AxisAlignedBox aabb;
    int i;
    bool found = false;

    for (i = 0; i < this->nodes.size(); ++i) {
        if (!(this->nodes[i].isParent)) {
            if (leafIdx == 1) {
                aabb = {
                    glm::vec3(this->nodes[i].bounds[0][0], this->nodes[i].bounds[1][0], this->nodes[i].bounds[2][0]), 
                    glm::vec3(this->nodes[i].bounds[0][1], this->nodes[i].bounds[1][1], this->nodes[i].bounds[2][1])
                };
                found = true;
                break;
            } else if (leafIdx <= 0)
                break;
            else
                --leafIdx;
        }
    }

    if (found) {
        // Draw the AABB as a transparent green box.
        // AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
        // drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

        // Draw the AABB as a (white) wireframe box.
        // drawAABB(aabb, DrawMode::Wireframe);
        drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

        // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
        for (int j = 0; j < this->nodes[i].indexes.size(); ++++j) {
            glm::uvec3 triangle = m_pScene->meshes[this->nodes[i].indexes[j]].triangles[this->nodes[i].indexes[j+1]];
            drawTriangle(m_pScene->meshes[this->nodes[i].indexes[j]].vertices[triangle.x], 
                m_pScene->meshes[this->nodes[i].indexes[j]].vertices[triangle.y], 
                m_pScene->meshes[this->nodes[i].indexes[j]].vertices[triangle.z]);
        }
    }
}


// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    
    
    Ray nZero;
    Ray nOne ;
    Ray nTwo;
    Ray nPoint ;
    Vertex hv0; // vertices of hit triangle
    Vertex hv1; // needed to draw visual debug for normal interpolation and traversal
    Vertex hv2;

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
                    glm::vec3 barc = computeBarycentricCoord(v0.position, v1.position, v2.position, ray.origin + ray.t * ray.direction);
                    glm::vec2 texC = barc.x * v0.texCoord + barc.y * v1.texCoord + barc.z * v2.texCoord;
                    hitInfo.barycentricCoord = barc;
                    hitInfo.texCoord = texC;
                    
                    if (features.enableNormalInterp) { //interpolate normal and update information to draw normals for visual debug
                        
                        
                        glm::vec3 interpolatedNormal = interpolateNormal(v0.normal, v1.normal, v2.normal, barc);
                        hitInfo.normal = normalize(interpolatedNormal);
                        nZero = Ray { v0.position, v0.normal, 1 };
                        nOne = Ray { v1.position, v1.normal, 1 };
                        nTwo = Ray { v2.position, v2.normal, 1 };
                        nPoint = Ray { ray.origin + ray.t * ray.direction, hitInfo.normal, 1 };
                        
                        
                        
                    } else {
                        glm::vec3 d1 = v1.position - v0.position;
                        glm::vec3 d2 = v2.position - v0.position;

                        glm::vec3 normal = normalize(glm::cross(d1, d2));
                        hitInfo.normal = normal;
                    }
                    
                    
                }
            }
        }
       
        
        if (hit == true && features.enableNormalInterp) { // if a hit occurs and normal interpolation enabled, draw normals
            drawRay(nZero, glm::vec3(1, 0, 0));
            drawRay(nOne, glm::vec3(1, 0, 0));
            drawRay(nTwo, glm::vec3(1, 0, 0));
            drawRay(nPoint, glm::vec3(0, 0, 1));
        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        bool hit = false;

        std::priority_queue<pqNode, std::vector<pqNode>, pqNode> pq; // priority queue sorts for the smalles t value of pqNode

        pqNode current(nodes[0]);

        pq.push(current);

        while (!pq.empty()) {

            current = pq.top();
            pq.pop();
            glm::vec3 lower(current.node.bounds[0][0], current.node.bounds[1][0], current.node.bounds[2][0]);
            glm::vec3 upper(current.node.bounds[0][1], current.node.bounds[1][1], current.node.bounds[2][1]);
            AxisAlignedBox box(lower, upper);

            drawAABB(box, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
            if (!current.node.isParent) {

                for (int i = 0; i < current.node.indexes.size(); i = i + 2) {
                    Mesh mesh = m_pScene->meshes[current.node.indexes[i]];

                    const auto triangle = mesh.triangles[current.node.indexes[i + 1]];
                    Vertex v0 = mesh.vertices[triangle[0]];
                    Vertex v1 = mesh.vertices[triangle[1]];
                    Vertex v2 = mesh.vertices[triangle[2]];

                    if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                        hitInfo.material = mesh.material;
                        hit = true;
                        glm::vec3 barc = computeBarycentricCoord(v0.position, v1.position, v2.position, ray.origin + ray.t * ray.direction);
                        glm::vec2 texC = barc.x * v0.texCoord + barc.y * v1.texCoord + barc.z * v2.texCoord;
                        hitInfo.barycentricCoord = barc;
                        hitInfo.texCoord = texC;
                        hv0 = v0;
                        hv1 = v1;
                        hv2 = v2;
                        if (features.enableNormalInterp) { // interpolate normal and update information to draw normals for visual debug

                            glm::vec3 interpolatedNormal = interpolateNormal(v0.normal, v1.normal, v2.normal, barc);
                            hitInfo.normal = normalize(interpolatedNormal);
                            nZero = Ray { v0.position, v0.normal, 1 };
                            nOne = Ray { v1.position, v1.normal, 1 };
                            nTwo = Ray { v2.position, v2.normal, 1 };
                            nPoint = Ray { ray.origin + ray.t * ray.direction, hitInfo.normal, 1 };

                        } else {
                            glm::vec3 d1 = v1.position - v0.position;
                            glm::vec3 d2 = v2.position - v0.position;

                            glm::vec3 normal = normalize(glm::cross(d1, d2));
                            hitInfo.normal = normal;
                        }
                    }
                }

                if (pq.empty() || ray.t <= pq.top().t) {
                    break;
                }

            } else {
                for (int i : current.node.indexes) {
                    pqNode child = pqNode(nodes[i]);
                    if (intersectRayPQNode(ray, child)) {
                        pq.push(child);
                    }
                }
            }
        }

        if (!pq.empty() && features.enableShading) { // draw the nodes which were already enqued but not visited when shading and bvh are enabled
            while (!pq.empty()) {
                current = pq.top();
                pq.pop();
                glm::vec3 lower(current.node.bounds[0][0], current.node.bounds[1][0], current.node.bounds[2][0]);
                glm::vec3 upper(current.node.bounds[0][1], current.node.bounds[1][1], current.node.bounds[2][1]);
                AxisAlignedBox box(lower, upper);

                drawAABB(box, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 0.05f), 0.1f);
            }
        }
        if (hit) {
            drawTriangle(hv0, hv1, hv2);
            if (features.enableNormalInterp) { // if a hit occurs and normal interpolation enabled, draw normals
                drawRay(nZero, glm::vec3(1, 0, 0));
                drawRay(nOne, glm::vec3(1, 0, 0));
                drawRay(nTwo, glm::vec3(1, 0, 0));
                drawRay(nPoint, glm::vec3(0, 0, 1));
            }
            return true;
        }
        return hit;
    }
}

// returns if pqNode intersect with Ray
// only returns true if t of node >0 and t of node < t of ray
bool BoundingVolumeHierarchy::intersectRayPQNode(Ray& ray, pqNode& node) const
{

    glm::vec3 lower(node.node.bounds[0][0], node.node.bounds[1][0], node.node.bounds[2][0]);
    glm::vec3 upper(node.node.bounds[0][1], node.node.bounds[1][1], node.node.bounds[2][1]);
    AxisAlignedBox box(lower, upper);
    float tmax = FLT_MAX;
    float tmin = FLT_MIN;
    if (ray.direction.x != 0) {
        float tx_min = (box.lower.x - ray.origin.x) / ray.direction.x;
        float tx_max = (box.upper.x - ray.origin.x) / ray.direction.x;
        tmin = fmax(tmin, fmin(tx_min, tx_max));
        tmax = fmin(tmax, fmax(tx_min, tx_max));
    } else {
        if (ray.origin.x < box.lower.x || ray.origin.x > box.upper.x) {
            return false;
        }
    }

    if (ray.direction.x != 0) {
        float ty_min = (box.lower.y - ray.origin.y) / ray.direction.y;
        float ty_max = (box.upper.y - ray.origin.y) / ray.direction.y;
        tmin = fmax(tmin, fmin(ty_min, ty_max));
        tmax = fmin(tmax, fmax(ty_min, ty_max));
    } else {
        if (ray.origin.y < box.lower.y || ray.origin.y > box.upper.y) {
            return false;
        }
    }

    if (ray.direction.z != 0) {
        float tz_min = (box.lower.z - ray.origin.z) / ray.direction.z;
        float tz_max = (box.upper.z - ray.origin.z) / ray.direction.z;
        tmin = fmax(tmin, fmin(tz_min, tz_max));
        tmax = fmin(tmax, fmax(tz_min, tz_max));
    } else {
        if (ray.origin.z < box.lower.z || ray.origin.z > box.upper.z) {
            return false;
        }
    }
    if (tmax < tmin) {
        return false;
    }
    if (tmin > 0 && ray.t > tmin) {
        node.t = tmin;
        return true;
    }
    return false;
}
