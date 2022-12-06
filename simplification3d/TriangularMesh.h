#ifndef _TRIANGULAR_MESH_H_
#define _TRIANGULAR_MESH_H_

#include <cwgUtilities.h>
#include "Cleaverlib/Cleaver.h"
#include "TspVideo.h"
#include "Tetrahedron.h"
#include "TetEdgeHeap.h"
#include "TetCube.h"

namespace cwg
{
	class TriangularMesh
	{
	public:
		TriangularMesh() { TetTriangle::reset_nextid(); }
		~TriangularMesh();

		inline TetTriangle*							GetTri(int k) const { return Tris[k]; }
		bool										save_ply(const std::string& fn, TetrahedralMesh* tetmesh) const;

		std::vector<TetTriangle*>					Tris;
	private:

	};
}

#endif