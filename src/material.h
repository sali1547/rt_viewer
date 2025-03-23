#pragma once
#ifndef MATERIALH
#define MATERIALH
#include <glm/glm.hpp>

struct HitRecord;

#include "rt_ray.h"
#include "rt_hitable.h"






namespace rt{
    glm::vec3 reflect(const glm::vec3& v, const glm::vec3& n) {
        return v -2*glm::dot(v,n)*n;
    }

    glm::vec3 random_in_unit_sphere() {
		glm::vec3 p;
		do {
			p = (float)2.0*glm::vec3(random(), random(), random()) - glm::vec3(1, 1, 1);
		} while ((p.x*p.x + p.y*p.y + p.z*p.z) >= 1.0);
		return p;
	}

    class material {
        public:
            virtual bool scatter(const Ray& r_in, const HitRecord& rec, glm::vec3& attenuation, Ray& scattered) const = 0;
    }; 


    class lambertian : public material {
        public:
            lambertian(const glm::vec3& a) : albedo(a) {}
            virtual bool scatter (const Ray& r_in, const HitRecord& rec, glm::vec3& attenuation, Ray& scattered) const {
                glm::vec3 target = rec.p + rec.normal + random_in_unit_sphere();
                scattered = Ray(rec.p, target-rec.p);
                attenuation = albedo; 
                return true;
            }
            glm::vec3 albedo;
    };

class metal : public material {
    public:
        metal(const glm::vec3& a, float f) : albedo(a), fuzz(glm::clamp(f, 0.0f, 1.0f)) {}
        virtual bool scatter(const Ray& r_in, const HitRecord& rec, glm::vec3& attenuation, Ray& scattered) 
        const override {
            glm::vec3 reflected = reflect(glm::normalize(r_in.direction()), rec.normal);
            scattered = Ray(rec.p, reflected + fuzz * random_in_unit_sphere()); 
            attenuation = albedo;
            return (glm::dot(scattered.direction(), rec.normal) > 0); 
        }
        glm::vec3 albedo;
        float fuzz;
        };
}
#endif 


