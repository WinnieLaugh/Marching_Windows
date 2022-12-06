#ifndef _VEC4_H_
#define _VEC4_H_

#include <iostream>
#include "../SymMat4.h"
namespace Cleaver
{
	class vec4
	{
	public:
		vec4(double x=0, double y=0, double z=0, double g=0) : x(x), y(y), z(z), g(g) {}

	public:
		double						x;
		double						y;
		double						z;
		double						g;

		vec4&						operator=(const vec4 &a);
		vec4&						operator+=(const vec4 &a);
		vec4&						operator*=(double c);
		vec4&						operator/=(double c);

		double&						operator[](const size_t);
		double						operator[](const size_t) const;

		inline double				dot(const vec4 &b) const { return this->x*b.x + this->y*b.y + this->z*b.z + this->g*b.g; }
		inline void					normalize() { double len = sqrt(x*x + y*y + z*z + g*g)+std::numeric_limits<double>::epsilon();
													x /= len; y /= len; z /= len; g /= len; }

		static vec4					zero;
		static vec4					unitX;
		static vec4					unitY;
		static vec4					unitZ;
		static vec4					unitG;

		static vec4					min(const vec4 &a, const vec4 &b);
		static vec4					max(const vec4 &a, const vec4 &b);

		cwg::SymMat4&				vvt();

		std::string					toString() const;

		friend std::ostream&		operator<<(std::ostream &stream, const vec4 &v);
	};

	double							dot(const vec4 &a, const vec4 &b);
	double							length(const vec4 &a);
	double							L1(const vec4 &a);
	double							L2(const vec4 &a);
	vec4							normalize(const vec4 &v1);

	vec4							operator+(const vec4 &a, const vec4 &b);
	vec4							operator-(const vec4 &a, const vec4 &b);
	vec4							operator*(const vec4 &a, double b);
	vec4							operator*(double a, const vec4 &b);
	vec4							operator/(const vec4 &a, double b);
}

#endif