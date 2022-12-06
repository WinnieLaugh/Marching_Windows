#include "stdafx.h"
#include "vec4.h"
#include "math.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

namespace Cleaver
{
	// static variables
	vec4 vec4::zero(0,0,0);
	vec4 vec4::unitX(1,0,0);
	vec4 vec4::unitY(0,1,0);
	vec4 vec4::unitZ(0,0,1);

	vec4& vec4::operator=( const vec4 &a )
	{
		this->x = a.x;
		this->y = a.y;
		this->z = a.z;
		this->g = a.g;

		return *this;
	}

	vec4& vec4::operator+=( const vec4 &a )
	{
		this->x += a.x;
		this->y += a.y;
		this->z += a.z;
		this->g += a.g;

		return *this;
	}

	vec4& vec4::operator*=( double c )
	{
		this->x *= c;
		this->y *= c;
		this->z *= c;
		this->g *= c;

		return *this;
	}

	vec4& vec4::operator/=( double c )
	{
		this->x /= c;
		this->y /= c;
		this->z /= c;
		this->g /= c;

		return *this;
	}

	double& vec4::operator[]( const size_t idx)
	{
		switch(idx)
		{
		case 0: return this->x;
		case 1: return this->y;
		case 2: return this->z;
		case 3: return this->g;
		default: throw -1;  // bad index
		}
	}

	double vec4::operator[]( const size_t idx) const
	{
		switch(idx)
		{
		case 0: return this->x;
		case 1: return this->y;
		case 2: return this->z;
		case 3: return this->g;
		default: throw -1;  // bad index
		}
	}

	Cleaver::vec4 vec4::min( const vec4 &a, const vec4 &b )
	{
		return vec4((a.x < b.x) ? a.x : b.x,
					(a.y < b.y) ? a.y : b.y,
					(a.z < b.z) ? a.z : b.z,
					(a.g < b.g) ? a.g : b.g);
	}

	Cleaver::vec4 vec4::max( const vec4 &a, const vec4 &b )
	{
		return vec4((a.x > b.x) ? a.x : b.x,
					(a.y > b.y) ? a.y : b.y,
					(a.z > b.z) ? a.z : b.z,
					(a.g > b.g) ? a.g : b.g);
	}

	std::string vec4::toString() const
	{
		std::stringstream ss;
		ss << "[" << std::setprecision(5) << this->x << ", " << this->y << ", " << this->z << ", " << this->g << "]";
		return ss.str();
	}

	std::ostream &operator<<(std::ostream &stream, const vec4 &v)
	{
		stream << std::fixed;
		return stream << std::setprecision(3) << v.x << " " << v.y << " " << v.z << " " << v.g;
	}

	double dot( const vec4 &a, const vec4 &b )
	{
		return a.x*b.x + a.y*b.y + a.z*b.z + a.g*b.g;
	}

	double length( const vec4 &a )
	{
		return sqrt(a.x*a.x + a.y*a.y + a.z*a.z + a.g*a.g);
	}

	double L1( const vec4 &a )
	{
		return a.x + a.y + a.z + a.g;
	}

	double L2( const vec4 &a )
	{
		return sqrt(a.x*a.x + a.y*a.y + a.z*a.z + a.g*a.g);
	}

	Cleaver::vec4 normalize( const vec4 &v1 )
	{
		return v1 / length(v1);
	}

	Cleaver::vec4 operator+( const vec4 &a, const vec4 &b )
	{
		return vec4(a.x+b.x, a.y+b.y, a.z+b.z, a.g+b.g);
	}

	Cleaver::vec4 operator-( const vec4 &a, const vec4 &b )
	{
		return vec4(a.x-b.x, a.y-b.y, a.z-b.z, a.g-b.g);
	}

	Cleaver::vec4 operator*( const vec4 &a, double b )
	{
		return vec4(b*a.x, b*a.y, b*a.z, b*a.g);
	}

	Cleaver::vec4 operator*( double a, const vec4 &b )
	{
		return vec4(a*b.x, a*b.y, a*b.z, a*b.g);
	}

	Cleaver::vec4 operator/( const vec4 &a, double b )
	{
		return vec4(a.x/b, a.y/b, a.z/b, a.g/b);
	}

	cwg::SymMat4& vec4::vvt()
	{
		cwg::SymMat4 m;
		m.set2zero();

		double v[4] = {x,y,z,g};

		for (int i = 0; i < 4; i++)
		{
			for (int j=0; j<4; j++)
			{
				m(i,j) = v[i] * v[j];
			}
		}
		return m;
	}
}