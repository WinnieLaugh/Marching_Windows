#include "stdafx.h"
#include <cstdlib>
#include <cmath>
#include <fstream>
#include "BCCLattice3D.h"
#include "GeneralizedStencilTable.h"
#include "TetMesh.h"
#include "Volume.h"
#include "ScalarField.h"
#include "BCCLattice3DLite.h"

using namespace std;
using namespace Cleaver;

#pragma region extern_define
extern const int dualLongEdgeFaceGroup[6][4];
extern const int primalLongEdgeFaceGroup[12][4][2];
extern const int shortEdgeFaceGroup[8][6];
extern const int dualLongEdgeTetGroup[8][6];
extern const int primalLongEdgeTetGroup[12][4][2];
extern const int shortEdgeTetGroup[8][6];
extern const int faceVertexGroup[36][3][2];
extern const int faceEdgeGroup[36][3][2];
extern const int faceTetGroup[36][2];
extern const int edgeCellGroup[27][3];
extern const int vertexCellGroup[8][8][3];
extern const int vertexEdgeGroup[14][2];
extern const int vertexFaceGroup[36][2];
extern const int primalVertexFaceGroupVertices[36][6];
extern const int dualVertexFaceGroupVertices[36][6];
extern const int vertexTetGroup[24][2];
extern const int TetTetGroup[24][2];
extern const int FaceFaceGroup[36][2];
extern const int EdgeEdgeGroup[26][4][2];
extern const int VertexVertexGroup[8][8][4];

#ifndef CELL_ID
#define CELL_ID 0
#endif

#ifndef VERT_ID
#define VERT_ID 1
#endif

#ifndef EDGE_ID
#define EDGE_ID 1
#endif

#ifndef FACE_ID
#define FACE_ID 1
#endif

#ifndef TET_ID
#define TET_ID  1
#endif
#pragma endregion extern_define

Cleaver::BCCLattice3DLite::BCCLattice3DLite( const LabelVolume* tspvideo ) : tspvideo(tspvideo)
{
	m_iNumMaterials = tspvideo->nTsp();

	if (tspvideo->nPaddings == 2)
	{
		m_iWidth = tspvideo->Width();
		m_iHeight = tspvideo->Height();
		m_iDepth = tspvideo->nFrames();
	}
	else if (tspvideo->nPaddings == 4)
	{
		m_iWidth = (tspvideo->Width() + 1) / 2;
		m_iHeight = (tspvideo->Height() + 1) / 2;
		m_iDepth = (tspvideo->nFrames() + 1) / 2;
	}

	verts = new std::vector<Vertex3D*>();
	tets = new std::vector<Tet*>();

	this->tree = new Octree(m_iWidth-1,m_iHeight-1,m_iDepth-1, *verts, *tets);

	m_bDataLoaded = false;
	m_bCutsComputed = false;
	m_bTriplesComputed = false;
	m_bQuadsComputed = false;
	m_bGeneralized = false;
	m_bStenciled = false;
}


Cleaver::BCCLattice3DLite::~BCCLattice3DLite()
{
	// clean up any data on leaves of tree
	// must do this before we destroy the tree
	unsigned int total_cells = cut_cells.size() + buffer_cells.size();


	// must delete ALL tets, followed by ALL faces, followed by ALL edges, followed by ALL verts
	// reason behind this is generalization may have edges pointing to vertices, and if vertices
	// are deleted, edge->cut pointer is now invalid and could do goofy things.


	// delete tets
	for(unsigned int c=0; c < total_cells; c++)
	{
		OTCell *cell = NULL;
		if(c < cut_cells.size())
			cell = cut_cells[c];
		else
			cell = buffer_cells[c - cut_cells.size()];

		if(cell->tets)
		{
			for(int t=0; t < 24; t++)
			{
				// if tet exists
				if(cell->tets[t])
				{
					// set to NULL in adjacent cells that shares it
					int cell_index = TetTetGroup[t][CELL_ID];
					int tet_index  = TetTetGroup[t][TET_ID];
					OTCell *nCell = tree->getNeighbor(cell, edgeCellGroup[cell_index]);

					if(nCell && nCell->tets)
						nCell->tets[tet_index] = NULL;

					// delete tet
					delete cell->tets[t];
					cell->tets[t] = NULL;
				}
			}
			delete[] cell->tets;
			cell->tets = NULL;
		}
	}

	// delete faces
	for(unsigned int c=0; c < total_cells; c++)
	{
		OTCell *cell = NULL;
		if(c < cut_cells.size())
			cell = cut_cells[c];
		else
			cell = buffer_cells[c - cut_cells.size()];

		// delete faces
		if(cell->face)
		{
			for(int f=0; f < 36; f++)
			{
				// if face exists
				if(cell->face[f])
				{
					// first 12 faces are interior, just delete them, otherwise
					// set to NULL in adjacent cells that shares it first
					if(f >= 12)
					{
						int cell_index = FaceFaceGroup[f][CELL_ID];
						int face_index = FaceFaceGroup[f][FACE_ID];
						OTCell *nCell = tree->getNeighbor(cell, edgeCellGroup[cell_index]);

						if(nCell && nCell->face)
							nCell->face[face_index] = NULL;
					}

					// delete face
					delete cell->face[f];
					cell->face[f] = NULL;
				}
			}
			delete[] cell->face;
			cell->face = NULL;
		}
	}

	// delete edges
	for(unsigned int c=0; c < total_cells; c++)
	{
		OTCell *cell = NULL;
		if(c < cut_cells.size())
			cell = cut_cells[c];
		else
			cell = buffer_cells[c - cut_cells.size()];

		// delete edges
		if(cell->edge)
		{
			for(int e=0; e < 26; e++)
			{
				// if edge exists
				if(cell->edge[e])
				{
					Edge3D *edge = cell->edge[e];

					// set to NULL in adjacent cells that shares it
					if(e >= 8)
					{
						int group_count = e < 14 ? 2 : 4;

						for(int n=0; n < group_count; n++)
						{
							int cell_index = EdgeEdgeGroup[e][n][CELL_ID];
							int edge_index = EdgeEdgeGroup[e][n][EDGE_ID];
							OTCell *nCell = tree->getNeighbor(cell, edgeCellGroup[cell_index]);

							if(nCell && nCell->edge)
								nCell->edge[edge_index] = NULL;
						}
					}

					// delete edge
					delete edge;
					cell->edge[e] = NULL;
				}
			}

			delete[] cell->edge;
			cell->edge = NULL;
		}
	}


	// delete vertex 3ds
	for(unsigned int c=0; c < total_cells; c++)
	{
		OTCell *cell = NULL;
		if(c < cut_cells.size())
			cell = cut_cells[c];
		else
			cell = buffer_cells[c - cut_cells.size()];

		// delete vertex arrays
		if(cell->vert)
		{
			// Delete Primal Vertices
			for(int v=0; v < 8; v++)
			{
				if(cell->vert[v])
				{
					// save pointer to memory
					Vertex3D *vertex = cell->vert[v];

					// null out all 8 possible references
					for(int n=0; n < 8; n++)
					{
						const int *cell_offset = VertexVertexGroup[v][n];
						int vert_index = VertexVertexGroup[v][n][3];
						OTCell *nCell = tree->getNeighbor(cell, cell_offset);

						if(nCell && nCell->vert)
							nCell->vert[vert_index] = NULL;
					}

					// free memory if not used in output mesh
					if(vertex->tm_v_index < 0)
						delete vertex;
					cell->vert[v] = NULL;  // redundant safety
				}
			}


			// Delete Dual Vertex if not used in output Mesh
			if(cell->vert[C] && cell->vert[C]->tm_v_index < 0)
			{
				delete cell->vert[C];
				cell->vert[C] = NULL;
			}

			// Delete the Array
			delete[] cell->vert;
			cell->vert = NULL;
		}
	}

	// finally delete tree
	if(tree != NULL){
		delete tree;
		tree = NULL;
	}
}

int Cleaver::BCCLattice3DLite::Label3D( int iw, int jh, int kf ) const
{
	
	return ( (tspvideo->nPaddings == 2) ? (*tspvideo)(kf, jh, iw) : (*tspvideo)(kf*2, jh*2, iw*2) );
}

int Cleaver::BCCLattice3DLite::Label3D_C( int iw, int jh, int kf ) const
{
	return ( (tspvideo->nPaddings == 2) ? (*tspvideo)(kf, jh, iw) : (*tspvideo)(kf*2+1, jh*2+1, iw*2+1) );
}

void Cleaver::BCCLattice3DLite::getAdjacencyLists( const Face3D *face, Vertex3D *verts[3], Edge3D *edges[3] )
{
	getVertsAroundFace(face, verts);
	getEdgesAroundFace(face, edges);
}

void Cleaver::BCCLattice3DLite::getAdjacencyLists( const Tet3D *tet, Vertex3D *verts[4], Edge3D *edges[6], Face3D *faces[4] )
{
	OTCell *cell = tet->cell;

	switch(tet->tet_index)
	{
		// right tets
	case TRU:
		{
			// upper
			OTCell *Rcell = tree->getNeighbor(cell, 1, 0, 0);

			verts[0] = cell->vert[C];
			verts[1] = cell->vert[URB];
			verts[2] = Rcell->vert[C];
			verts[3] = cell->vert[URF];

			edges[0] = cell->edge[DURB];
			edges[1] = cell->edge[CR];
			edges[2] = cell->edge[URF];
			edges[3] = Rcell->edge[DULB];
			edges[4] = Rcell->edge[DULF];
			edges[5] = cell->edge[UR];

			faces[0] = cell->face[FRUB];
			faces[1] = cell->face[FRUF];
			faces[2] = cell->face[FUR];
			faces[3] = Rcell->face[FUL];

			break;
		}
	case TRL:
		{
			// lower
			OTCell *Rcell = tree->getNeighbor(cell, 1, 0, 0);

			verts[0] = cell->vert[C];
			verts[1] = cell->vert[LRF];
			verts[2] = Rcell->vert[C];
			verts[3] = cell->vert[LRB];

			edges[0] = cell->edge[DLRF];
			edges[1] = cell->edge[CR];
			edges[2] = cell->edge[DLRB];
			edges[3] = Rcell->edge[DLLF];
			edges[4] = Rcell->edge[DLLB];
			edges[5] = cell->edge[LR];

			faces[0] = cell->face[FRLF];
			faces[1] = cell->face[FRLB];
			faces[2] = cell->face[FLR];
			faces[3] = Rcell->face[FLL];

			break;
		}
	case TRF:
		{
			// front
			OTCell *Rcell = tree->getNeighbor(cell, 1, 0, 0);

			verts[0] = cell->vert[C];
			verts[1] = cell->vert[URF];
			verts[2] = Rcell->vert[C];
			verts[3] = cell->vert[LRF];

			edges[0] = cell->edge[DURF];
			edges[1] = cell->edge[CR];
			edges[2] = cell->edge[DLRF];
			edges[3] = Rcell->edge[DULF];
			edges[4] = Rcell->edge[DLLF];
			edges[5] = cell->edge[FR];

			faces[0] = cell->face[FRUF];
			faces[1] = cell->face[FRLF];
			faces[2] = cell->face[FFR];
			faces[3] = Rcell->face[FFL];

			break;
		}
	case TRB:
		{
			// back
			OTCell *Rcell = tree->getNeighbor(cell, 1, 0, 0);

			verts[0] = cell->vert[C];
			verts[1] = cell->vert[LRB];
			verts[2] = Rcell->vert[C];
			verts[3] = cell->vert[URB];

			edges[0] = cell->edge[DLRB];
			edges[1] = cell->edge[CR];
			edges[2] = cell->edge[DURB];
			edges[3] = Rcell->edge[DLLB];
			edges[4] = Rcell->edge[DULB];
			edges[5] = cell->edge[BR];

			faces[0] = cell->face[FRLB];
			faces[1] = cell->face[FRUB];
			faces[2] = cell->face[FBR];
			faces[3] = Rcell->face[FBL];

			break;
		}
		// upper tets
	case TUF:
		{
			// front
			OTCell *Ucell = tree->getNeighbor(cell, 0, 1, 0);

			verts[0] = cell->vert[C];
			verts[1] = cell->vert[ULF];
			verts[2] = Ucell->vert[C];
			verts[3] = cell->vert[URF];

			edges[0] = cell->edge[DULF];
			edges[1] = cell->edge[CU];
			edges[2] = cell->edge[DURF];
			edges[3] = Ucell->edge[DLLF];
			edges[4] = Ucell->edge[DLRF];
			edges[5] = cell->edge[UF];

			faces[0] = cell->face[FUFL];
			faces[1] = cell->face[FUFR];
			faces[2] = cell->face[FUF];
			faces[3] = Ucell->face[FLF];

			break;
		}
	case TUB:
		{
			// back
			OTCell *Ucell = tree->getNeighbor(cell, 0, 1, 0);

			verts[0] = cell->vert[C];
			verts[1] = cell->vert[URB];
			verts[2] = Ucell->vert[C];
			verts[3] = cell->vert[ULB];

			edges[0] = cell->edge[DURB];
			edges[1] = cell->edge[CU];
			edges[2] = cell->edge[DULB];
			edges[3] = Ucell->edge[DLRB];
			edges[4] = Ucell->edge[DLLB];
			edges[5] = cell->edge[UB];

			faces[0] = cell->face[FUBR];
			faces[1] = cell->face[FUBL];
			faces[2] = cell->face[FUB];
			faces[3] = Ucell->face[FLB];

			break;
		}
	case TUL:
		{
			// left
			OTCell *Ucell = tree->getNeighbor(cell, 0, 1, 0);

			verts[0] = cell->vert[C];
			verts[1] = cell->vert[ULB];
			verts[2] = Ucell->vert[C];
			verts[3] = cell->vert[ULF];

			edges[0] = cell->edge[DULB];
			edges[1] = cell->edge[CU];
			edges[2] = cell->edge[DULF];
			edges[3] = Ucell->edge[DLLB];
			edges[4] = Ucell->edge[DLLF];
			edges[5] = cell->edge[UL];

			faces[0] = cell->face[FUBL];
			faces[1] = cell->face[FUFL];
			faces[2] = cell->face[FUL];
			faces[3] = Ucell->face[FLL];

			break;
		}
	case TUR:
		{
			// right
			OTCell *Ucell = tree->getNeighbor(cell, 0, 1, 0);

			verts[0] = cell->vert[C];
			verts[1] = cell->vert[URF];
			verts[2] = Ucell->vert[C];
			verts[3] = cell->vert[URB];

			edges[0] = cell->edge[DURF];
			edges[1] = cell->edge[CU];
			edges[2] = cell->edge[DURB];
			edges[3] = Ucell->edge[DLRF];
			edges[4] = Ucell->edge[DLRB];
			edges[5] = cell->edge[UR];

			faces[0] = cell->face[FUFR];
			faces[1] = cell->face[FUBR];
			faces[2] = cell->face[FUR];
			faces[3] = Ucell->face[FLR];

			break;
		}
		// back tets
	case TBT:
		{
			// top
			OTCell * Bcell = tree->getNeighbor(cell, 0, 0, 1);

			verts[0] = cell->vert[C];
			verts[1] = cell->vert[ULB];
			verts[2] = Bcell->vert[C];
			verts[3] = cell->vert[URB];

			edges[0] = cell->edge[DULB];
			edges[1] = cell->edge[CB];
			edges[2] = cell->edge[DURB];
			edges[3] = Bcell->edge[DULF];
			edges[4] = Bcell->edge[DURF];
			edges[5] = cell->edge[UB];

			faces[0] = cell->face[FBUL];
			faces[1] = cell->face[FBUR];
			faces[2] = cell->face[FUB];
			faces[3] = Bcell->face[FUF];

			break;
		}
	case TBB:
		{
			// bottom
			OTCell * Bcell = tree->getNeighbor(cell, 0, 0, 1);

			verts[0] = cell->vert[C];
			verts[1] = cell->vert[LRB];
			verts[2] = Bcell->vert[C];
			verts[3] = cell->vert[LLB];

			edges[0] = cell->edge[DLRB];
			edges[1] = cell->edge[CB];
			edges[2] = cell->edge[DLLB];
			edges[3] = Bcell->edge[DLRF];
			edges[4] = Bcell->edge[DLLF];
			edges[5] = cell->edge[LB];

			faces[0] = cell->face[FBLR];
			faces[1] = cell->face[FBLL];
			faces[2] = cell->face[FLB];
			faces[3] = Bcell->face[FLF];

			break;
		}
	case TBL:
		{
			// left
			OTCell * Bcell = tree->getNeighbor(cell, 0, 0, 1);

			verts[0] = cell->vert[C];
			verts[1] = cell->vert[LLB];
			verts[2] = Bcell->vert[C];
			verts[3] = cell->vert[ULB];

			edges[0] = cell->edge[DLLB];
			edges[1] = cell->edge[CB];
			edges[2] = cell->edge[DULB];
			edges[3] = Bcell->edge[DLLF];
			edges[4] = Bcell->edge[DULF];
			edges[5] = cell->edge[BL];

			faces[0] = cell->face[FBLL];
			faces[1] = cell->face[FBUL];
			faces[2] = cell->face[FBL];
			faces[3] = Bcell->face[FFL];

			break;
		}
	case TBR:
		{
			// right
			OTCell * Bcell = tree->getNeighbor(cell, 0, 0, 1);

			verts[0] = cell->vert[C];
			verts[1] = cell->vert[URB];
			verts[2] = Bcell->vert[C];
			verts[3] = cell->vert[LRB];

			edges[0] = cell->edge[DURB];
			edges[1] = cell->edge[CB];
			edges[2] = cell->edge[DLRB];
			edges[3] = Bcell->edge[DURF];
			edges[4] = Bcell->edge[DLRF];
			edges[5] = cell->edge[BR];

			faces[0] = cell->face[FBUR];
			faces[1] = cell->face[FBLR];
			faces[2] = cell->face[FBR];
			faces[3] = Bcell->face[FFR];

			break;
		}
	}
}


void Cleaver::BCCLattice3DLite::getVertsAroundFace( const Face3D *face, Vertex3D *verts[3] )
{
	for(int i=0; i < 3; i++)
	{
		int cell_id = faceVertexGroup[face->face_index][i][CELL_ID];
		int vert_id = faceVertexGroup[face->face_index][i][VERT_ID];
		OTCell *cell = tree->getNeighbor(face->cell, edgeCellGroup[cell_id]);
		verts[i] = cell->vert[vert_id];

	}
}

void Cleaver::BCCLattice3DLite::getEdgesAroundFace( const Face3D *face, Edge3D *edges[3] )
{
	for(int i=0; i < 3; i++)
	{
		int cell_id = faceEdgeGroup[face->face_index][i][CELL_ID];
		int edge_id = faceEdgeGroup[face->face_index][i][EDGE_ID];
		OTCell *cell = tree->getNeighbor(face->cell, edgeCellGroup[cell_id]);
		edges[i] = cell->edge[edge_id];
	}
}

unsigned char Cleaver::BCCLattice3DLite::keyFromAdjacentEdges( Edge3D *edges[6] )
{
	unsigned char key = 0; // 64;
	(edges[0]->cut && edges[0]->cut->order() == CUT) ? key |= 32 : 0;
	(edges[1]->cut && edges[1]->cut->order() == CUT) ? key |= 16 : 0 ;
	(edges[2]->cut && edges[2]->cut->order() == CUT) ? key |=  8 : 0 ;
	(edges[3]->cut && edges[3]->cut->order() == CUT) ? key |=  4 : 0 ;
	(edges[4]->cut && edges[4]->cut->order() == CUT) ? key |=  2 : 0 ;
	(edges[5]->cut && edges[5]->cut->order() == CUT) ? key |=  1 : 0 ;

	return key;
}

bool Cleaver::BCCLattice3DLite::isKeyValid( unsigned char key )
{
	if (key ==  0 || key == 11 || key == 22 || key == 29 || key == 31 ||
		key == 37 || key == 46 || key == 47 ||
		key == 51 || key == 55 || key == 56 || key == 59 || key == 61 || key == 62 ||
		key == 63)
		return true;
	else
		return false;
}

void Cleaver::BCCLattice3DLite::getRightHandedVertexList( const Tet3D *tet, Vertex3D *verts[15] )
{
	OTCell *cell = tet->cell;

	switch(tet->tet_index)
	{
		// right tets
	case TRU:
		{
			// upper
			OTCell *Rcell = tree->getNeighbor(cell, 1, 0, 0);

			verts[   _A] = cell->vert[C];
			verts[   _B] = cell->vert[URB];
			verts[   _C] = Rcell->vert[C];
			verts[   _D] = cell->vert[URF];
			verts[  _AB] = cell->edge[DURB]->cut;
			verts[  _AC] = cell->edge[CR]->cut;
			verts[  _AD] = cell->edge[URF]->cut;
			verts[  _CD] = Rcell->edge[DULF]->cut;
			verts[  _BD] = cell->edge[UR]->cut;
			verts[  _BC] = Rcell->edge[DULB]->cut;
			verts[ _BCD] = Rcell->face[FUL]->triple;
			verts[ _ACD] = cell->face[FRUF]->triple;
			verts[ _ABD] = cell->face[FUR]->triple;
			verts[ _ABC] = cell->face[FRUB]->triple;
			verts[_ABCD] = tet->quad;
			break;
		}
	case TRL:
		{
			// lower
			OTCell *Rcell = tree->getNeighbor(cell, 1, 0, 0);

			verts[   _A] = cell->vert[C];
			verts[   _B] = cell->vert[LRF];
			verts[   _C] = Rcell->vert[C];
			verts[   _D] = cell->vert[LRB];
			verts[  _AB] = cell->edge[DLRF]->cut;
			verts[  _AC] = cell->edge[CR]->cut;
			verts[  _AD] = cell->edge[DLRB]->cut;
			verts[  _CD] = Rcell->edge[DLLB]->cut;
			verts[  _BD] = cell->edge[LR]->cut;
			verts[  _BC] = Rcell->edge[DLLF]->cut;
			verts[ _BCD] = Rcell->face[FLL]->triple;
			verts[ _ACD] = cell->face[FRLB]->triple;
			verts[ _ABD] = cell->face[FLR]->triple;
			verts[ _ABC] = cell->face[FRLF]->triple;
			verts[_ABCD] = tet->quad;
			break;
		}
	case TRF:
		{
			// front
			OTCell *Rcell = tree->getNeighbor(cell, 1, 0, 0);

			verts[   _A] = cell->vert[C];
			verts[   _B] = cell->vert[URF];
			verts[   _C] = Rcell->vert[C];
			verts[   _D] = cell->vert[LRF];
			verts[  _AB] = cell->edge[DURF]->cut;
			verts[  _AC] = cell->edge[CR]->cut;
			verts[  _AD] = cell->edge[DLRF]->cut;
			verts[  _CD] = Rcell->edge[DLLF]->cut;
			verts[  _BD] = cell->edge[FR]->cut;
			verts[  _BC] = Rcell->edge[DULF]->cut;
			verts[ _BCD] = Rcell->face[FFL]->triple;
			verts[ _ACD] = cell->face[FRLF]->triple;
			verts[ _ABD] = cell->face[FFR]->triple;
			verts[ _ABC] = cell->face[FRUF]->triple;
			verts[_ABCD] = tet->quad;
			break;
		}
	case TRB:
		{
			// back
			OTCell *Rcell = tree->getNeighbor(cell, 1, 0, 0);

			verts[   _A] = cell->vert[C];
			verts[   _B] = cell->vert[LRB];
			verts[   _C] = Rcell->vert[C];
			verts[   _D] = cell->vert[URB];
			verts[  _AB] = cell->edge[DLRB]->cut;
			verts[  _AC] = cell->edge[CR]->cut;
			verts[  _AD] = cell->edge[DURB]->cut;
			verts[  _CD] = Rcell->edge[DULB]->cut;
			verts[  _BD] = cell->edge[BR]->cut;
			verts[  _BC] = Rcell->edge[DLLB]->cut;
			verts[ _BCD] = Rcell->face[FBL]->triple;
			verts[ _ACD] = cell->face[FRUB]->triple;
			verts[ _ABD] = cell->face[FBR]->triple;
			verts[ _ABC] = cell->face[FRLB]->triple;
			verts[_ABCD] = tet->quad;
			break;
		}
		// upper tets
	case TUF:
		{
			// front
			OTCell *Ucell = tree->getNeighbor(cell, 0, 1, 0);

			verts[   _A] = cell->vert[C];
			verts[   _B] = cell->vert[ULF];
			verts[   _C] = Ucell->vert[C];
			verts[   _D] = cell->vert[URF];
			verts[  _AB] = cell->edge[DULF]->cut;
			verts[  _AC] = cell->edge[CU]->cut;
			verts[  _AD] = cell->edge[DURF]->cut;
			verts[  _CD] = Ucell->edge[DLRF]->cut;
			verts[  _BD] = cell->edge[UF]->cut;
			verts[  _BC] = Ucell->edge[DLLF]->cut;
			verts[ _BCD] = Ucell->face[FLF]->triple;
			verts[ _ACD] = cell->face[FUFR]->triple;
			verts[ _ABD] = cell->face[FUF]->triple;
			verts[ _ABC] = cell->face[FUFL]->triple;
			verts[_ABCD] = tet->quad;
			break;
		}
	case TUB:
		{
			// back
			OTCell *Ucell = tree->getNeighbor(cell, 0, 1, 0);

			verts[   _A] = cell->vert[C];
			verts[   _B] = cell->vert[URB];
			verts[   _C] = Ucell->vert[C];
			verts[   _D] = cell->vert[ULB];
			verts[  _AB] = cell->edge[DURB]->cut;
			verts[  _AC] = cell->edge[CU]->cut;
			verts[  _AD] = cell->edge[DULB]->cut;
			verts[  _CD] = Ucell->edge[DLLB]->cut;
			verts[  _BD] = cell->edge[UB]->cut;
			verts[  _BC] = Ucell->edge[DLRB]->cut;
			verts[ _BCD] = Ucell->face[FLB]->triple;
			verts[ _ACD] = cell->face[FUBL]->triple;
			verts[ _ABD] = cell->face[FUB]->triple;
			verts[ _ABC] = cell->face[FUBR]->triple;
			verts[_ABCD] = tet->quad;
			break;
		}
	case TUL:
		{
			// left
			OTCell *Ucell = tree->getNeighbor(cell, 0, 1, 0);
			verts[   _A] = cell->vert[C];
			verts[   _B] = cell->vert[ULB];
			verts[   _C] = Ucell->vert[C];
			verts[   _D] = cell->vert[ULF];
			verts[  _AB] = cell->edge[DULB]->cut;
			verts[  _AC] = cell->edge[CU]->cut;
			verts[  _AD] = cell->edge[DULF]->cut;
			verts[  _CD] = Ucell->edge[DLLF]->cut;
			verts[  _BD] = cell->edge[UL]->cut;
			verts[  _BC] = Ucell->edge[DLLB]->cut;
			verts[ _BCD] = Ucell->face[FLL]->triple;
			verts[ _ACD] = cell->face[FUFL]->triple;
			verts[ _ABD] = cell->face[FUL]->triple;
			verts[ _ABC] = cell->face[FUBL]->triple;
			verts[_ABCD] = tet->quad;
			break;
		}
	case TUR:
		{
			// right
			OTCell *Ucell = tree->getNeighbor(cell, 0, 1, 0);

			verts[   _A] = cell->vert[C];
			verts[   _B] = cell->vert[URF];
			verts[   _C] = Ucell->vert[C];
			verts[   _D] = cell->vert[URB];
			verts[  _AB] = cell->edge[DURF]->cut;
			verts[  _AC] = cell->edge[CU]->cut;
			verts[  _AD] = cell->edge[DURB]->cut;
			verts[  _CD] = Ucell->edge[DLRB]->cut;
			verts[  _BD] = cell->edge[UR]->cut;
			verts[  _BC] = Ucell->edge[DLRF]->cut;
			verts[ _BCD] = Ucell->face[FLR]->triple;
			verts[ _ACD] = cell->face[FUBR]->triple;
			verts[ _ABD] = cell->face[FUR]->triple;
			verts[ _ABC] = cell->face[FUFR]->triple;
			verts[_ABCD] = tet->quad;
			break;
		}
		// back tets
	case TBT:
		{
			// top
			OTCell * Bcell = tree->getNeighbor(cell, 0, 0, 1);

			verts[   _A] = cell->vert[C];
			verts[   _B] = cell->vert[ULB];
			verts[   _C] = Bcell->vert[C];
			verts[   _D] = cell->vert[URB];
			verts[  _AB] = cell->edge[DULB]->cut;
			verts[  _AC] = cell->edge[CB]->cut;
			verts[  _AD] = cell->edge[DURB]->cut;
			verts[  _CD] = Bcell->edge[DURF]->cut;
			verts[  _BD] = cell->edge[UB]->cut;
			verts[  _BC] = Bcell->edge[DULF]->cut;
			verts[ _BCD] = Bcell->face[FUF]->triple;
			verts[ _ACD] = cell->face[FBUR]->triple;
			verts[ _ABD] = cell->face[FUB]->triple;
			verts[ _ABC] = cell->face[FBUL]->triple;
			verts[_ABCD] = tet->quad;
			break;
		}
	case TBB:
		{
			// bottom
			OTCell * Bcell = tree->getNeighbor(cell, 0, 0, 1);

			verts[   _A] = cell->vert[C];
			verts[   _B] = cell->vert[LRB];
			verts[   _C] = Bcell->vert[C];
			verts[   _D] = cell->vert[LLB];
			verts[  _AB] = cell->edge[DLRB]->cut;
			verts[  _AC] = cell->edge[CB]->cut;
			verts[  _AD] = cell->edge[DLLB]->cut;
			verts[  _CD] = Bcell->edge[DLLF]->cut;
			verts[  _BD] = cell->edge[LB]->cut;
			verts[  _BC] = Bcell->edge[DLRF]->cut;
			verts[ _BCD] = Bcell->face[FLF]->triple;
			verts[ _ACD] = cell->face[FBLL]->triple;
			verts[ _ABD] = cell->face[FLB]->triple;
			verts[ _ABC] = cell->face[FBLR]->triple;
			verts[_ABCD] = tet->quad;
			break;
		}
	case TBL:
		{
			// left
			OTCell * Bcell = tree->getNeighbor(cell, 0, 0, 1);

			verts[   _A] = cell->vert[C];
			verts[   _B] = cell->vert[LLB];
			verts[   _C] = Bcell->vert[C];
			verts[   _D] = cell->vert[ULB];
			verts[  _AB] = cell->edge[DLLB]->cut;
			verts[  _AC] = cell->edge[CB]->cut;
			verts[  _AD] = cell->edge[DULB]->cut;
			verts[  _CD] = Bcell->edge[DULF]->cut;
			verts[  _BD] = cell->edge[BL]->cut;
			verts[  _BC] = Bcell->edge[DLLF]->cut;
			verts[ _BCD] = Bcell->face[FFL]->triple;
			verts[ _ACD] = cell->face[FBUL]->triple;
			verts[ _ABD] = cell->face[FBL]->triple;
			verts[ _ABC] = cell->face[FBLL]->triple;
			verts[_ABCD] = tet->quad;
			break;
		}
	case TBR:
		{
			// right
			OTCell * Bcell = tree->getNeighbor(cell, 0, 0, 1);

			verts[   _A] = cell->vert[C];
			verts[   _B] = cell->vert[URB];
			verts[   _C] = Bcell->vert[C];
			verts[   _D] = cell->vert[LRB];
			verts[  _AB] = cell->edge[DURB]->cut;
			verts[  _AC] = cell->edge[CB]->cut;
			verts[  _AD] = cell->edge[DLRB]->cut;
			verts[  _CD] = Bcell->edge[DLRF]->cut;
			verts[  _BD] = cell->edge[BR]->cut;
			verts[  _BC] = Bcell->edge[DURF]->cut;
			verts[ _BCD] = Bcell->face[FFR]->triple;
			verts[ _ACD] = cell->face[FBLR]->triple;
			verts[ _ABD] = cell->face[FBR]->triple;
			verts[ _ABC] = cell->face[FBUR]->triple;
			verts[_ABCD] = tet->quad;
			break;
		}
	default:
		{
			cerr << "Fatal Error: InvalidTet Index!!" << endl;
			exit(51);
		}
	}
}

void Cleaver::BCCLattice3DLite::getEdgesAroundTet( const Tet3D *tet, Edge3D *edges[6] )
{
	OTCell *cell = tet->cell;

	switch(tet->tet_index)
	{
		// right tets
	case TRU:
		{
			// upper
			OTCell *Rcell = tree->getNeighbor(cell, 1, 0, 0);

			edges[0] = cell->edge[DURB];
			edges[1] = cell->edge[CR];
			edges[2] = cell->edge[URF];
			edges[3] = Rcell->edge[DULB];
			edges[4] = Rcell->edge[DULF];
			edges[5] = cell->edge[UR];

			break;
		}
	case TRL:
		{
			// lower
			OTCell *Rcell = tree->getNeighbor(cell, 1, 0, 0);

			edges[0] = cell->edge[DLRF];
			edges[1] = cell->edge[CR];
			edges[2] = cell->edge[DLRB];
			edges[3] = Rcell->edge[DLLF];
			edges[4] = Rcell->edge[DLLB];
			edges[5] = cell->edge[LR];

			break;
		}
	case TRF:
		{
			// front
			OTCell *Rcell = tree->getNeighbor(cell, 1, 0, 0);

			edges[0] = cell->edge[DURF];
			edges[1] = cell->edge[CR];
			edges[2] = cell->edge[DLRF];
			edges[3] = Rcell->edge[DULF];
			edges[4] = Rcell->edge[DLLF];
			edges[5] = cell->edge[FR];

			break;
		}
	case TRB:
		{
			// back
			OTCell *Rcell = tree->getNeighbor(cell, 1, 0, 0);

			edges[0] = cell->edge[DLRB];
			edges[1] = cell->edge[CR];
			edges[2] = cell->edge[DURB];
			edges[3] = Rcell->edge[DLLB];
			edges[4] = Rcell->edge[DULB];
			edges[5] = cell->edge[BR];

			break;
		}
		// upper tets
	case TUF:
		{
			// front
			OTCell *Ucell = tree->getNeighbor(cell, 0, 1, 0);

			edges[0] = cell->edge[DULF];
			edges[1] = cell->edge[CU];
			edges[2] = cell->edge[DURF];
			edges[3] = Ucell->edge[DLLF];
			edges[4] = Ucell->edge[DLRF];
			edges[5] = cell->edge[UF];

			break;
		}
	case TUB:
		{
			// back
			OTCell *Ucell = tree->getNeighbor(cell, 0, 1, 0);

			edges[0] = cell->edge[DURB];
			edges[1] = cell->edge[CU];
			edges[2] = cell->edge[DULB];
			edges[3] = Ucell->edge[DLRB];
			edges[4] = Ucell->edge[DLLB];
			edges[5] = cell->edge[UB];

			break;
		}
	case TUL:
		{
			// left
			OTCell *Ucell = tree->getNeighbor(cell, 0, 1, 0);

			edges[0] = cell->edge[DULB];
			edges[1] = cell->edge[CU];
			edges[2] = cell->edge[DULF];
			edges[3] = Ucell->edge[DLLB];
			edges[4] = Ucell->edge[DLLF];
			edges[5] = cell->edge[UL];

			break;
		}
	case TUR:
		{
			// right
			OTCell *Ucell = tree->getNeighbor(cell, 0, 1, 0);

			edges[0] = cell->edge[DURF];
			edges[1] = cell->edge[CU];
			edges[2] = cell->edge[DURB];
			edges[3] = Ucell->edge[DLRF];
			edges[4] = Ucell->edge[DLRB];
			edges[5] = cell->edge[UR];

			break;
		}
		// back tets
	case TBT:
		{
			// top
			OTCell * Bcell = tree->getNeighbor(cell, 0, 0, 1);

			edges[0] = cell->edge[DULB];
			edges[1] = cell->edge[CB];
			edges[2] = cell->edge[DURB];
			edges[3] = Bcell->edge[DULF];
			edges[4] = Bcell->edge[DURF];
			edges[5] = cell->edge[UB];

			break;
		}
	case TBB:
		{
			// bottom
			OTCell * Bcell = tree->getNeighbor(cell, 0, 0, 1);

			edges[0] = cell->edge[DLRB];
			edges[1] = cell->edge[CB];
			edges[2] = cell->edge[DLLB];
			edges[3] = Bcell->edge[DLRF];
			edges[4] = Bcell->edge[DLLF];
			edges[5] = cell->edge[LB];

			break;
		}
	case TBL:
		{
			// left
			OTCell * Bcell = tree->getNeighbor(cell, 0, 0, 1);

			edges[0] = cell->edge[DLLB];
			edges[1] = cell->edge[CB];
			edges[2] = cell->edge[DULB];
			edges[3] = Bcell->edge[DLLF];
			edges[4] = Bcell->edge[DULF];
			edges[5] = cell->edge[BL];

			break;
		}
	case TBR:
		{
			// right
			OTCell * Bcell = tree->getNeighbor(cell, 0, 0, 1);

			edges[0] = cell->edge[DURB];
			edges[1] = cell->edge[CB];
			edges[2] = cell->edge[DLRB];
			edges[3] = Bcell->edge[DURF];
			edges[4] = Bcell->edge[DLRF];
			edges[5] = cell->edge[BR];

			break;
		}
	}
}
