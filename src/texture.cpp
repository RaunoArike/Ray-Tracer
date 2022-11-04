#include "texture.h"
#include <framework/image.h>
#include <glm/common.hpp>
#include <iostream>
#include <cmath>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    
    if (features.extra.enableBilinearTextureFiltering) {
        float x = texCoord.x * image.width - 0.5f;
        float y = (1 - texCoord.y) * image.height - 0.5f;
        int iu0 = std::floor(x);
        int iu1 = iu0 + 1;
        int iv0 = std::floor(y);
        int iv1 = iv0 + 1;
        float a_u = iu1 - x;
        float b_u = 1 - a_u;
        float a_v = iv1 - y;
        float b_v = 1 - a_v;
        iu0 = glm::clamp(iu0, 0, image.width - 1);
        iu1 = glm::clamp(iu1, 0, image.width - 1);
        iv0 = glm::clamp(iv0, 0, image.height - 1);
        iv1 = glm::clamp(iv1, 0, image.height - 1);
        return a_u * a_v * image.pixels[iu0 * image.width + iv0] + a_u * b_v * image.pixels[iu0 * image.width + iv1] + b_u * a_v * image.pixels[iu1 * image.width + iv0] + b_u * b_v * image.pixels[iu1 * image.width + iv1];
    } else {
        int x = (texCoord.x * image.width);
        int y = (texCoord.y * image.height);
        int i = (image.height - y - 1) * image.width + x;
        return image.pixels[i];
    }
}


