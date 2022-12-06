#ifndef _CWG_BCCLATTICE3D_H_
#define _CWG_BCCLATTICE3D_H_

#include "basic_geometry.h"

#include <cwgUtilities.h>
#include "../TSPVideo.h"

namespace Cleaver
{
	class AbstractVolume;
	class ScalarField;
	class Poly3D;
	class Tet3D;
	class Face3D;
	class Edge3D;
	class Octree;
	class OTCell;
	class BCCLattice3DLite;
	class Tet;

	class BCCLattice3DLite
	{
	public:
		BCCLattice3DLite(const LabelVolume* tspvideo);
		~BCCLattice3DLite();

		const LabelVolume*									tspvideo;
		Octree*												tree;		  // octree structure

		std::vector<OTCell *>								cut_cells;    // cells that actually have cuts in them
		std::vector<OTCell *>								buffer_cells; // cells that are adjacent to ones with cuts

		int													Label3D(int iw, int jh, int kf) const;
		int													Label3D_C(int iw, int jh, int kf) const;

		inline int											materials() const { return m_iNumMaterials; }
		inline int											width()     const { return m_iWidth; }
		inline int											height()    const { return m_iHeight; }
		inline int											depth()     const { return m_iDepth; }
		
		unsigned char										keyFromAdjacentEdges(Edge3D *edges[6]);
		bool												isKeyValid(unsigned char key);
		
		void												getRightHandedVertexList(const Tet3D *tet, Vertex3D *verts[15]);

		void												getVertsAroundFace(const Face3D *face, Vertex3D *verts[3]);
		void												getEdgesAroundFace(const Face3D *face, Edge3D *edges[3]);
		
		void												getEdgesAroundTet(const Tet3D *tet, Edge3D *edges[6]);

		void												getAdjacencyLists(const Face3D *face, Vertex3D *verts[3], Edge3D *edges[3]);
		void												getAdjacencyLists(const Tet3D *tet, Vertex3D *verts[4], Edge3D *edges[6], Face3D *faces[4]);

		std::vector<Vertex3D*>*								verts;
		std::vector<Tet*>*									tets;

		//-------------------------------------------- State Setters ------------------------------------------------//
		void												setDataLoaded(bool state)		{ m_bDataLoaded = state; }
		void												setCutsComputed(bool state)		{ m_bCutsComputed = state; }
		void												setTriplesComputed(bool state)	{ m_bTriplesComputed = state; }
		void												setQuadsComputed(bool state)	{ m_bQuadsComputed = state; }
		void												setGeneralized(bool state)		{ m_bGeneralized = state; }
		void												setStenciled(bool state)		{ m_bStenciled = state; }

	private:
		int													m_iNumMaterials;    // number of materials in volume
		int													m_iWidth;           // width of grid
		int													m_iHeight;          // height of grid
		int													m_iDepth;           // depth of grid

		//-------------------------------------------- State of Data -----------------------------------------------//
		bool												m_bDataLoaded;
		bool												m_bCutsComputed;
		bool												m_bTriplesComputed;
		bool												m_bQuadsComputed;
		bool												m_bGeneralized;
		bool												m_bStenciled;
	};


}

#endif