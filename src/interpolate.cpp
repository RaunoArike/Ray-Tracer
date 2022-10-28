#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord (const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    
    glm::vec3 n = normalize(cross(v1 - v0, v2 - v0));
    float ABC = dot(n, cross(v1 - v0, v2 - v0));
    float PBC = dot(n, cross(v1 - p, v2 - p));
    float PCA = dot(n, cross(v2 - p, v0 - p));
    
    float alpha = PBC / ABC;
    float beta = PCA / ABC;
    float gamma = 1 - alpha - beta;
   
    
    return glm::vec3(alpha,beta,gamma);
}

glm::vec3 interpolateNormal (const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
     
    glm::vec3 interNormal = barycentricCoord.x * n0 + barycentricCoord.y * n1 + barycentricCoord.z * n2;
    return interNormal;
    
}

glm::vec2 interpolateTexCoord (const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
// TODO: implement this function.
    return barycentricCoord.x * t0 + barycentricCoord.y * t1 + barycentricCoord.z * t2  ;
}
