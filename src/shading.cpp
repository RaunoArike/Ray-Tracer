#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{

    glm::vec3 light = normalize((lightPosition - (ray.origin + ray.t * ray.direction)));
    glm::vec3 view = normalize(-ray.direction);
    
    glm::vec3 res(0.0);
    if (dot(hitInfo.normal, light) < 0) {
        return glm::vec3(0, 0, 0);
    }

    if (dot(hitInfo.normal, view) < 0) {
        return glm::vec3(0, 0, 0);
    }
    if (features.enableTextureMapping &&hitInfo.material.kdTexture) {
        Image tex = *hitInfo.material.kdTexture;
        res = lightColor * (acquireTexel(tex,hitInfo.texCoord,features) * dot(normalize(hitInfo.normal), light));
    } else {
        res = lightColor * (hitInfo.material.kd * dot(normalize(hitInfo.normal), light));
    }
        

    return res;
    
  
    
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    Ray reflectionRay {};
    hitInfo.normal = glm::normalize(hitInfo.normal);
    reflectionRay.direction = glm::normalize(ray.direction - 2 * glm::dot(ray.direction, hitInfo.normal) * hitInfo.normal);
    auto offset = 0.001f * hitInfo.normal;
    reflectionRay.origin = ray.origin + ray.t * ray.direction + offset;
    return reflectionRay;
}