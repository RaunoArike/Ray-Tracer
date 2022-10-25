#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord (const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    glm::vec3 a = v1 - v0;
    glm::vec3 b = v2 - v0;
    glm::vec3 c = p - v0;
    float aa = dot(a, a);
    float ab = dot(a, b);
    float bb = dot(b, b);
    float ca = dot(c, a);
    float cb = dot(c, a);
    float d = aa * bb - ab * ab;
    float v = (bb * ca - ab * cb) / d;
    float w = (aa * cb - ab * ca) / d;
    float u = 1.0f - v - w;
    
    return glm::vec3(u,v,w);
}

glm::vec3 interpolateNormal (const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
     
    glm::vec3 interNormal = barycentricCoord.x * n0 + barycentricCoord.y * n1 + barycentricCoord.z * n2;
    return interNormal;
    
}

glm::vec2 interpolateTexCoord (const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
// TODO: implement this function.
    return glm::vec2(0.0);
}
