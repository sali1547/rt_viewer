#pragma once

#include "rt_hitable.h"

namespace rt {

class Sphere : public Hitable {
  public:
    Sphere() {}
    Sphere(const glm::vec3 &cen, float r) : center(cen), radius(r) {};
	Sphere(const glm::vec3 &cen, float r, material* mat) : center(cen), radius(r) , mat_ptr(mat) {};
    virtual bool hit(const Ray &r, float t_min, float t_max, HitRecord &rec) const;
    glm::vec3 center;
    float radius;
    material *mat_ptr;
};

// Ray-sphere test from "Ray Tracing in a Weekend" book (page 16)
bool Sphere::hit(const Ray &r, float t_min, float t_max, HitRecord &rec) const {
    glm::vec3 oc = r.origin() - center;
    float a = glm::dot(r.direction(), r.direction());
    float b = glm::dot(oc, r.direction());
    float c = glm::dot(oc, oc) - radius * radius;
    float discriminant = b * b - a * c;

    if (discriminant > 0) {
        float temp = (-b - sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr;  // Store material in hit record
            return true;
        }
        temp = (-b + sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr;  // Store material in hit record
            return true;
        }
    }
    return false;
}

}  // namespace rt
