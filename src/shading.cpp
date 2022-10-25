#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    
    
    glm::vec3 light = normalize((lightPosition-(ray.origin+ray.t*ray.direction)));
    if (dot(hitInfo.normal, light) < 0) {
        return glm::vec3(0, 0, 0);
    }
    glm::vec3 view = normalize(-ray.direction);
    if (dot(hitInfo.normal, view) < 0) {
        return glm::vec3(0, 0, 0);
    }
    
    glm::vec3 reflct = normalize(glm::reflect(-light, normalize(hitInfo.normal)));
    
    return lightColor * hitInfo.material.kd * dot(normalize(hitInfo.normal), light) + lightColor * hitInfo.material.ks * pow(dot(reflct, view), hitInfo.material.shininess);
  
    
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    Ray reflectionRay {};
    // TODO: implement the reflection ray computation.
    return reflectionRay;
}