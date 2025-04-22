#pragma once
#include "types.h"

class Plane {
public:
    float d;
    Vector3f n;

    Plane(const Vector3f& _nd) : d(_nd.norm()), n(_nd / d) {}
    Plane(const Vector3f& _n, float _d) : d(_d), n(_n) {}
    Plane(const Vector3f& a, const Vector3f& b, const Vector3f& c) {
        n = (b - a).cross(c - a).normalized();
        d = a.dot(n);
    }

    bool intersect(const Vector3f& a, const Vector3f& b, Vector3f& result) const {
        const float ad = a.dot(n) - d;
        const float bd = b.dot(n) - d;
        if (ad * bd >= 0.0f) { return false; }
        result = (ad * b - bd * a) / (ad - bd);
        return true;
    }

    float signed_distance(Vector3f p) const {
        return p.dot(n) - d;
    }
};
using Planes = std::vector<Plane>;
