#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>
#include <random>
#include <texture.cpp>


// returns a random number from a uniform distribution
float random(float range_start, float range_end)
{
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<float> distribution(range_start, range_end);
    return distribution(generator);
}

// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    float a = random(0.0f, 1.0f);
    position = a * segmentLight.endpoint0 + (1-a) * segmentLight.endpoint1;
    color = a * segmentLight.color0 + (1-a) * segmentLight.color1;
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    float a = random(0.0f, 1.0f);
    float b = random(0.0f, 1.0f);
    position = parallelogramLight.v0 + a * parallelogramLight.edge01 + b * parallelogramLight.edge02;
    color = a * b * parallelogramLight.color3 + (1-a) * b * parallelogramLight.color2 + a * (1-b) * parallelogramLight.color1 + (1 - a) * (1 - b) * parallelogramLight.color0;
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    // Compute intersection point
    glm::vec3 intersection = ray.origin + ray.t * ray.direction;
    float t = glm::distance(intersection, samplePos);

    // Compute ray from intersection to light source
    Ray intersectionToLight { intersection + glm::normalize(samplePos - intersection) * 0.001f, glm::normalize(samplePos - intersection), t }; // add small offset factor to reduce noise

    // Use bvh to find any triangles between light source and intersection
    bool result = bvh.intersect(intersectionToLight, hitInfo, features);

    if (result) {
        drawRay(intersectionToLight, { 1, 0, 0 }); // Red indicates no direct light from this source
        return 0.0;
    } else {
        // reset origin and t for drawing the ray
        intersectionToLight.origin = intersection;
        intersectionToLight.t = t;
        drawRay(intersectionToLight, debugColor); // Passed colour indicates direct light
    }
    return 1.0;     // Return 1.0 if there is direct light from the sample position.
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 res = glm::vec3(0.0);
    glm::vec3 position { 0.0 };
    glm::vec3 color { 0.0 };
    if (features.enableShading) {
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                position = pointLight.position;
                color = pointLight.color;
                res += computeShading(position, color, features, ray, hitInfo);
                // if hard shadows are enabled, check whether the sample is visible
                if (features.enableHardShadow) {
                    res *= testVisibilityLightSample(position, color, bvh, features, ray, hitInfo);
                }
            } else if (std::holds_alternative<SegmentLight>(light)) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                glm::vec3 samplePos { 0.0f };
                glm::vec3 sampleColor { 0.0f };
                if (features.enableSoftShadow) {
                    float segment_sample_size = 40;
                    for (int i = 0; i < segment_sample_size; i++) {
                        // calculate sample position and color
                        sampleSegmentLight(segmentLight, samplePos, sampleColor);
                        // add the sample to the color if it's visible
                        res += computeShading(samplePos, sampleColor, features, ray, hitInfo) * testVisibilityLightSample(samplePos, sampleColor, bvh, features, ray, hitInfo);
                    }
                    // average the shading result from all rays that hit the light source
                    res /= segment_sample_size;
                }
            } else if (std::holds_alternative<ParallelogramLight>(light)) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                glm::vec3 samplePos { 0.0f };
                glm::vec3 sampleColor { 0.0f };
                if (features.enableSoftShadow) {
                    float paral_sample_size = 100;
                    for (int i = 0; i < paral_sample_size; i++) {
                        // calculate sample position and color
                        sampleParallelogramLight(parallelogramLight, samplePos, sampleColor);
                        // add the sample to the color if it's visible
                        res += computeShading(samplePos, sampleColor, features, ray, hitInfo) * testVisibilityLightSample(samplePos, sampleColor, bvh, features, ray, hitInfo);
                    }
                    // average the shading result from all rays that hit the light source
                    res /= paral_sample_size;
                }
            }
        }
        return res;

    } else {
        // If shading is disabled, return the albedo of the material.
        if (features.enableTextureMapping && hitInfo.material.kdTexture) {
            Image tex = *hitInfo.material.kdTexture;
            return acquireTexel(tex, hitInfo.texCoord, features) ;
        }
        return hitInfo.material.kd;
        
        
    }
}
