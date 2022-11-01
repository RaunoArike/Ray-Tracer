#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>


bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
    // TODO: implement this function.
    return false;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    // TODO: implement this function.
    return false;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    // TODO: implement this function.
    Plane plane;
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    // TODO: implement this function.
    return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    // TODO: implement this function.
    return false;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool getT(const AxisAlignedBox& box, Ray& ray)
{
    float tmax = FLT_MAX;
    float tmin = FLT_MIN;
    if (ray.direction.x != 0) {
        float tx_min = (box.lower.x - ray.origin.x) / ray.direction.x;
        float tx_max = (box.upper.x - ray.origin.x) / ray.direction.x;
        tmin = fmax(tmin, fmin(tx_min, tx_max));
        tmax = fmin(tmax, fmax(tx_min, tx_max));
    } else {
        if (ray.origin.x < box.lower.x || ray.origin.x > box.upper.x) {
            return FLT_MAX;
        }
    }

    if (ray.direction.x != 0) {
        float ty_min = (box.lower.y - ray.origin.y) / ray.direction.y;
        float ty_max = (box.upper.y - ray.origin.y) / ray.direction.y;
        tmin = fmax(tmin, fmin(ty_min, ty_max));
        tmax = fmin(tmax, fmax(ty_min, ty_max));
    } else {
        if (ray.origin.y < box.lower.y || ray.origin.y > box.upper.y) {
            return FLT_MAX;
        }
    }

    if (ray.direction.z != 0) {
        float tz_min = (box.lower.z - ray.origin.z) / ray.direction.z;
        float tz_max = (box.upper.z - ray.origin.z) / ray.direction.z;
        tmin = fmax(tmin, fmin(tz_min, tz_max));
        tmax = fmin(tmax, fmax(tz_min, tz_max));
    } else {
        if (ray.origin.z < box.lower.z || ray.origin.z > box.upper.z) {
            return FLT_MAX;
        }
    }
    if (tmax < tmin) {
        return FLT_MAX;
    }
    if (ray.t > tmin && tmin > 0) {
        
        return tmin;
    }
    
    return FLT_MAX;
    
}

// returns the
bool intersectAABB(AxisAlignedBox box, Ray ray)
{
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

    if (tmin > 0) {

        return true;
    }
    return false;
}


