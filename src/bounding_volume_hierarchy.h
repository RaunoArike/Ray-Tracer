#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>

// Forward declaration.
struct Scene;

class BoundingVolumeHierarchy {
private:
    struct Node {
        bool isParent; // Type; false implies a leaf node
        std::vector < std::vector<float>> bounds;    // { {x_min,x_max},{y_min,y_max},{z_min,z_max} } 
        std::vector<int> indexes; // Child node indexes XOR mesh + triangle indexes. [mesh0, triangle0, mesh1, triangle1, mesh2, triangle2, ...]
    };
    struct pqNode {
        Node node; //
        float t; // the t when the ray enters the box

        bool operator()(const pqNode& lhs, const pqNode& rhs) // Comperator used to sort the priority queue
        {
            return lhs.t > rhs.t; // gives smallest t the highest priority
        }
    };


public:
    // Helper function for constructor. Grows the tree using recursion.
    void growBVH(int nodeIndex, int recursionDepth);

    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level, Node node);
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;

    // Return true if ray hits the AABB of the priorityQueue Node
    // Only finds hits closer then the t stored
    // sets the t for the pqNode
    bool intersectRayPQNode(Ray& ray, pqNode& node) const;



private:
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
    std::vector<Node> nodes;                                    // Vector of nodes within the BVH.

    const int maxLevels = 10;                                   // Maximum depth of the BVH !!! 
};