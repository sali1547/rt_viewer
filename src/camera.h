#ifndef CAMERAH
#define CAMERAH

#include "rt_ray.h"


class camera {
    public:
        camera() {
            lower_left_corner = glm::vec3(-2.0, -1.0, -1.0);
            horizontal = glm::vec3(4.0, 0.0, 0.0);
            vertical = glm::vec3(0.0, 2.0, 0.0);
            origin = glm::vec3(0.0, 0.0, 0.0);
        }
        rt::Ray get_ray(float u, float v) { return rt::Ray(origin, lower_left_corner + u*horizontal + v*vertical - origin); }

        glm::vec3 origin;
        glm::vec3 lower_left_corner;
        glm::vec3 horizontal;
        glm::vec3 vertical; 
};

#endif
