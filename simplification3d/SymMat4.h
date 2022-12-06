#ifndef _SYMMAT7_H_
#define _SYMMAT7_H_

#include <cmath>
#include <iostream>

namespace cwg
{
	class SymMat4
	{
	public:
		SymMat4();
		~SymMat4() {}

		SymMat4			operator+( const SymMat4& B) const;
		SymMat4&		operator+=( const SymMat4& B);
		SymMat4			operator*( double val ) const;
		SymMat4&		operator*=( double val);

		double			operator()( int i, int j ) const;
		double&			operator()( int i, int j );

		inline void		set2zero() { memset(m_data, 0, sizeof(double)*16); }

	private:
		double			m_data[16];
	};

	std::ostream& operator<< (std::ostream& out, const SymMat4& symmat4);
}

#endif