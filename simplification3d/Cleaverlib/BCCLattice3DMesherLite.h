#ifndef _CWG_BCCLATTICE3D_MESHER_LITE_H_
#define _CWG_BCCLATTICE3D_MESHER_LITE_H_

#include "BCCLattice3DLite.h"
#include "TetMesh.h"

namespace Cleaver
{

class BCCLattice3DMesherLite
{
public:
	BCCLattice3DMesherLite(BCCLattice3DLite* lattice = NULL) : lattice(lattice) {}
	
	TetMesh*											mesh();

	void												compute_all_cuts();
	void												compute_all_trips();
	void												compute_all_quads();
	
	void												generalize_tets();
	void												fill_all_stencils();

	void												setLattice(BCCLattice3DLite* lattice){ this->lattice = lattice; }
	BCCLattice3DLite*									getLattice() { return this->lattice; }

protected:

private:

	BCCLattice3DLite*									lattice;
	TetMesh*											tm;

	void												compute_cut(Edge3D *edge);
	void												compute_triple(Face3D *face);
	void												compute_quadruple(Tet3D *tet);
	void												fill_stencil(Tet3D *tet);

	bool												isTransition(const Vertex3D* v1, const Vertex3D* v2);

	void												fixTriangleOrdering(Edge3D *edges[], Vertex3D *verts[]);
	void												fixTetrahedronOrdering(Face3D *faces[], Edge3D *edges[], Vertex3D *verts[]);
};

}

#endif