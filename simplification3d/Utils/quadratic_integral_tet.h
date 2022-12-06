#ifndef _QUADRATIC_INTEGRAL_TET_H_
#define _QUADRATIC_INTEGRAL_TET_H_

namespace hpan
{
	template<class T>
	class vec3g
	{
	public:
		vec3g<T>(T mx, T my, T mz) : x(mx), y(my), z(mz) {}
		T x;
		T y;
		T z;

		inline vec3g<T> operator- (const vec3g<T>& vm) const 
		{
			return vec3g<T> ( x - vm.x, y - vm.y, z - vm.z );
		}

		inline vec3g<T> operator+ (const vec3g<T>& vm) const
		{
			return vec3g<T> ( x + vm.x, y + vm.y, z + vm.z );
		}
	};

	template<class T> 
	T dot(vec3g<T> v1, vec3g<T> v2)
	{
		return ( v1.x * v2.x + v1.y * v2.y + v1.z * v2.z );
	}

	template <class T> inline T det3x3(
		T a00,  T a01,  T a02,
		T a10,  T a11,  T a12,
		T a20,  T a21,  T a22
		) 
	{
		T m01 = a00*a11 - a10*a01;
		T m02 = a00*a21 - a20*a01;
		T m12 = a10*a21 - a20*a11;
		return m01*a22 - m02*a12 + m12*a02;
	}

	template <class T> inline T mixed_product(
		const vec3g<T>& v1, const vec3g<T>& v2, const vec3g<T>& v3
		) 
	{
		return det3x3(
			v1.x, v2.x, v3.x,
			v1.y, v2.y, v3.y,
			v1.z, v2.z, v3.z
			) ;
	}

	template <class T> inline T tetra_signed_volume(
		const vec3g<T>& p, const vec3g<T>& q, 
		const vec3g<T>& r, const vec3g<T>& s
		) 
	{
		vec3g<T> U = q-p ;
		vec3g<T> V = r-p ;
		vec3g<T> W = s-p ;
		return mixed_product(U,V,W) / T(6) ;
	}

	template <class T> inline T quadratic_integral_tet(
		const vec3g<T>& seed, 
		const vec3g<T>& p1, const vec3g<T>& p2, const vec3g<T>& p3, const vec3g<T>& p4, 
		T& V
		) 
	{
		vec3g<T> v0 = p1 - seed;
		vec3g<T> v1 = p2 - seed;
		vec3g<T> v2 = p3 - seed;
		vec3g<T> v3 = p4 - seed;

		V = tetra_signed_volume(p1, p2, p3, p4);

		vec3g<T> v_23 = v2+v3;
		vec3g<T> v_123 = v1+v_23;
		vec3g<T> v_0123 = v0+v_123;
		return V * (dot(v0,v_0123)+dot(v1,v_123)+dot(v2,v_23)+dot(v3,v3)) / 10.;
	}
}

#endif
