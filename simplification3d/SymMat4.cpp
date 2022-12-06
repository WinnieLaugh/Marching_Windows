#include "stdafx.h"
#include "SymMat4.h"

cwg::SymMat4::SymMat4()
{
	memset(m_data, 0, sizeof(double)*16);
}

cwg::SymMat4 cwg::SymMat4::operator+( const SymMat4& B ) const
{
	SymMat4 rst;
	for (int i=0; i<16; i++)
		rst.m_data[i] = this->m_data[i] + B.m_data[i];

	return rst;
}

cwg::SymMat4& cwg::SymMat4::operator+=( const SymMat4& B )
{
	for (int i=0; i<16; i++)
		m_data[i] += B.m_data[i];

	return *this;
}

cwg::SymMat4 cwg::SymMat4::operator*( double val ) const
{
	SymMat4 rst;
	for (int i=0; i<16; i++)
		rst.m_data[i] = this->m_data[i] * val;
	return rst;
}

double cwg::SymMat4::operator()( int iRow, int iCol ) const
{
/*	if( iRow > iCol )
	{
		std::swap( iRow, iCol );
	}

	int iIndex = iRow*4 + iCol;
	for( int i = 0; i <= iRow; i++ )
	{
		iIndex -= i;
	}
	return m_data[iIndex];
	*/
	return m_data[iRow * 4 + iCol];

}

double& cwg::SymMat4::operator()( int iRow, int iCol )
{
/*	if( iRow > iCol )
	{
		std::swap( iRow, iCol );
	}

	int iIndex = iRow*4 + iCol;
	for( int i = 0; i <= iRow; i++ )
	{
		iIndex -= i;
	}
	return m_data[iIndex];
*/
	return m_data[iRow * 4 + iCol];
}

cwg::SymMat4& cwg::SymMat4::operator*=( double val )
{
	for (int i=0; i<16; i++)
		m_data[i] *= val;
	return *this;
}

std::ostream& cwg::operator<<( std::ostream& out, const SymMat4& symmat4 )
{
	for (int i=0; i<4; i++)
	{
		for (int j=0; j<4; j++)
		{
			out << symmat4(i,j) << " ";
		}
		out << std::endl;
	}
	return out;
}
