#ifndef MATERIALH
#define MATERIALH
#include "rt_ray.h"
#include <glm/glm.hpp>
#include "rt_hitable.h"

struct HitRecord;

namespace rt{
    glm::vec3 reflect(const glm::vec3& v, const glm::vec3& n) {
        return v - 2*glm::dot(v,n)*n;
    }


class material {
    public:
        virtual bool scatter(const rt::Ray& r_in, const HitRecord& rec, glm::vec3& attenuation, rt::Ray& scattered) const = 0;
        
};

class lambertian : public material {
    public:
        lambertian(const glm::vec3& a) : albedo(a) {}
        virtual bool scatter(const rt::Ray& r_in, const HitRecord& rec, glm::vec3& attenuation, rt::Ray& scattered) const {
            glm::vec3 target = rec.p + rec.normal + rt::random_in_unit_sphere();
            scattered = rt::Ray(rec.p, target-rec.p);
            attenuation = albedo;
            return true;
        }
   
        glm::vec3 albedo;
};

class metal : public material {
    public:
        metal(const glm::vec3& a, float f) : albedo(a), fuzz(glm::clamp(f, 0.0f, 1.0f)) {}

        virtual bool scatter(const rt::Ray& r_in, const HitRecord& rec, glm::vec3& attenuation, rt::Ray& scattered) const override {
            glm::vec3 reflected = rt::reflect(glm::normalize(r_in.direction()), rec.normal);
            scattered = rt::Ray(rec.p, reflected + fuzz * rt::random_in_unit_sphere()); // Add fuzziness
            attenuation = albedo;
            return (glm::dot(scattered.direction(), rec.normal) > 0);
        }

        glm::vec3 albedo;
        float fuzz;  // Controls roughness (0 = perfect mirror, 1 = rough)

};

}
#endif

