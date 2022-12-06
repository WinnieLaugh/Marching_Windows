#ifndef _TET_CUBE_H_
#define _TET_CUBE_H_

#include "TetVertex.h"
#include "Tetrahedron.h"

namespace cwg
{
	class TetCube
	{
	public:
		TetCube(int lb) : m_label(lb) {}
		~TetCube() {}

		int								Verts[8];

		void							ComputeTetrahedra(std::vector<Tetrahedron>::iterator& it);
		void							ComputeTetrahedra(std::vector<Tetrahedron>::iterator& it, int &tID);
		void							ComputeTetrahedra2(std::vector<Tetrahedron>::iterator& it, int &tID);
		void							ComputeTetrahedra3(std::vector<Tetrahedron>::iterator& it, int &tID);

	private:
		int								m_label;
	};
}

#endif