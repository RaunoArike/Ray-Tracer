#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>
#include <iostream>
#include <random>
#define EPSILON 0.0000001
#define MAX_DEPTH 20
#ifdef NDEBUG
#include <omp.h>
#endif

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {
        
        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.extra.enableTransparency) {
            // if alpha is less than one and ray depth is less than 20, do alpha blending
            if (hitInfo.material.transparency < 1 && rayDepth < MAX_DEPTH) {
                auto intersection = ray.origin + ray.direction * (ray.t + 0.0001f); // add a small offset to the ray
                auto rayContinuation = Ray { intersection, ray.direction };
                Lo = hitInfo.material.transparency * Lo + (1 - hitInfo.material.transparency) * getFinalColor(scene, bvh, rayContinuation, features, rayDepth + 1);
            }
        }

        if (features.enableRecursive) {
            Ray reflection = computeReflectionRay(ray, hitInfo);
            // TODO: put your own implementation of recursive ray tracing here.
            if (glm::length(hitInfo.material.ks) < EPSILON || rayDepth > MAX_DEPTH) {
                drawRay(ray, Lo);
                return Lo;
            }
            Lo += getFinalColor(scene, bvh, reflection, features, rayDepth+1);
        }

        if (features.enableNormalInterp) {

        }

        if (features.enableShading) {
            drawRay(ray, Lo);
        } else {
            // draw a white debug ray
            drawRay(ray, glm::vec3(1.0f));
        }

        // Set the color of the pixel to white if the ray hits.
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, depthOfFieldCalc(scene, bvh, cameraRay, features, camera, screen));
        }
    }
}

// returns a random number from a uniform distribution
float randNr(float range_start, float range_end)
{
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<float> distribution(range_start, range_end);
    return distribution(generator);
}

glm::vec3 depthOfFieldCalc(const Scene& scene, const BvhInterface& bvh, const Ray& ray, const Features& features, const Trackball& camera, Screen& screen)
{
    if (features.extra.enableDepthOfField) {
        auto focalLength = features.extra.enableDepthOfFieldF;
        auto aperture = features.extra.enableDepthOfFieldAperture;
        auto sampleCount = features.extra.enableDepthOfFieldSampleCount;
        
        auto focalPoint = ray.origin + focalLength * ray.direction;
        glm::vec3 color { 0.0f };

        for (int i = 0; i < sampleCount; i++) {
            // finding a random point on the aperture
            auto vecToPoint = glm::normalize(randNr(0.0f, 1.0f) * camera.up() + randNr(0.0f, 1.0f) * camera.left()) * aperture;
            vecToPoint.x = randNr(-1.0f, 1.0f) * vecToPoint.x;
            vecToPoint.y = randNr(-1.0f, 1.0f) * vecToPoint.y;
            vecToPoint.z = randNr(-1.0f, 1.0f) * vecToPoint.z;
            // vector from the origin of the camera ray to the randomly chosen point on the aperture
            auto startingPoint = ray.origin + vecToPoint;
            auto direction = glm::normalize(focalPoint - startingPoint);
            // ray from the point on the aperture to the focal point
            auto apertureRay = Ray { startingPoint, direction };
            color += getFinalColor(scene, bvh, apertureRay, features);
        }

        // average the result from all samples
        auto avgColor = color / float(sampleCount);

        // show the focal point for visual debug
        drawRay(Ray { ray.origin, ray.direction, focalLength }, avgColor);

        return avgColor;
    } else {
        return getFinalColor(scene, bvh, ray, features);
    }
}