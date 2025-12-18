#ifndef RTWEEKEND_H
#define RTWEEKEND_H
#include <cmath>
#include <limits>
#include <memory>
#include <random>

// Usings
using std::shared_ptr;
using std::make_shared;
using std::sqrt;

// Constants
const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// Utility Functions
inline double degrees_to_radians(double degrees) {
	return degrees * pi / 180.0;
}

inline double clamp(double x, double min, double max) {
	if (x < min) return min;
	if (x > max) return max;
	return x;
}

inline double random_double()
{
	static std::uniform_real_distribution<double> distribution(0.0, 1.0);
	static std::mt19937 generator;
	return distribution(generator);
}

inline double random_double(double min, double max)
{
	// Returns a random real in [min,max).
	return min + (max - min) * random_double();
}

inline vec3 random() {
	return vec3(random_double(), random_double(), random_double());
}

inline vec3 random(double min, double max) {
	return vec3(random_double(min, max), random_double(min, max),
		random_double(min, max));
}

vec3 random_in_unit_sphere() {
	while (true) {
		auto p = random(-1, 1);
		if (p.length_squared() >= 1) continue;
		return p;
	}
}

inline vec3 random_unit_vector() {
	return unit_vector(random_in_unit_sphere());
}

inline vec3 random_in_hemisphere(const vec3& normal) {
	vec3 in_unit_sphere = random_in_unit_sphere();
	if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
		return in_unit_sphere;
	else
		return -in_unit_sphere;
}

vec3 random_in_unit_disk() {
	while (true) {
		auto p = vec3(random_double(-1, 1), random_double(-1, 1), 0);
		if (p.length_squared() >= 1) continue;
		return p;
	}
}

// Common Headers
#include "ray.h"
#include "vec3.h"

#endif
