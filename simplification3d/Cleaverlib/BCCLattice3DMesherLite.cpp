#include "stdafx.h"
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include "BCCLattice3DMesherLite.h"
#include "GeneralizedInterfaceTable.h"
#include "GeneralizedMaterialTable.h"
#include "GeneralizedStencilTable.h"
#include "GeneralizedVertexTable.h"
#include "Volume.h"

using namespace std;
using namespace Cleaver;

extern int LATTICE_FACE_AXIS[36];

#ifndef CLEAVER_PI
#define CLEAVER_PI 3.14159265
#endif

TetMesh* Cleaver::BCCLattice3DMesherLite::mesh()
{
	compute_all_cuts();
	compute_all_trips();
	compute_all_quads();
	generalize_tets();
	
	tm = new TetMesh(*lattice->verts,*lattice->tets);
	fill_all_stencils();
	return tm;
}

void Cleaver::BCCLattice3DMesherLite::compute_all_cuts()
{
	for(unsigned int c=0; c < lattice->cut_cells.size(); c++)
	{
		OTCell *cell = lattice->cut_cells[c];

		for(int e=0; e < EDGES_PER_CELL; e++)
		{
			Edge3D *edge = cell->edge[e];
			if(!edge)
			{
				cerr << "Problem:  Material Transitions found on boundary." << endl;
				cerr << "Rerun with padding" << endl;
				exit(0);
			}

			if(!edge->evaluated)
				compute_cut(cell->edge[e]);
		}
	}

	lattice->setCutsComputed(true);
}

void Cleaver::BCCLattice3DMesherLite::compute_all_trips()
{
	for(unsigned int c = 0; c < lattice->cut_cells.size(); c++)
	{
		OTCell *cell = lattice->cut_cells[c];

		for(int f=0; f < FACES_PER_CELL; f++)
		{
			Face3D *face = cell->face[f];

			if(!face->evaluated)
				compute_triple(face);
		}
	}

	lattice->setTriplesComputed(true);
}

void Cleaver::BCCLattice3DMesherLite::compute_all_quads()
{
	for(unsigned int c = 0; c < lattice->cut_cells.size(); c++){

		OTCell *cell = lattice->cut_cells[c];

		for(int t=0; t < TETS_PER_CELL; t++)
		{
			Tet3D *tet = cell->tets[t];

			if(!tet->evaluated)
				compute_quadruple(tet);
		}
	}

	lattice->setQuadsComputed(true);
}

void Cleaver::BCCLattice3DMesherLite::generalize_tets()
{
	//--------------------------------------------------
	//      Generalize All Tets In Cut Cells
	//--------------------------------------------------
	for(unsigned int c = 0; c < lattice->cut_cells.size(); c++)
	{
		OTCell *cell = lattice->cut_cells[c];

		for(int t=0; t < TETS_PER_CELL; t++)
		{
			Tet3D *tet = cell->tets[t];

			// if no quad, start generalization
			if(tet && !tet->quad)
			{
				// look up generalization
				Vertex3D *verts[4];
				Edge3D *edges[6];
				Face3D *faces[4];

				lattice->getAdjacencyLists(tet, verts, edges, faces);

				unsigned char key = tet->key = lattice->keyFromAdjacentEdges(edges);

#ifdef  _DEBUG
				unsigned char key2 = lattice->keyFromAdjacentEdges(edges);
#endif

				if(!lattice->isKeyValid(key))
				{
#ifdef _DEBUG

					vec3 v0pos = verts[0]->pos();
					int	 v0lbl = verts[0]->label;
					vec3 v1pos = verts[1]->pos();
					int	 v1lbl = verts[1]->label;
					vec3 v2pos = verts[2]->pos();
					int	 v2lbl = verts[2]->label;
					vec3 v3pos = verts[3]->pos();
					int	 v3lbl = verts[3]->label;

					Edge3D* e0 = edges[0];
					vec3 e0v0 = e0->v1->pos();
					vec3 e0v1 = e0->v2->pos();
					Edge3D* e1 = edges[1];
					Edge3D* e2 = edges[2];
					Edge3D* e3 = edges[3];
					Edge3D* e4 = edges[4];
					Edge3D* e5 = edges[5];

					Face3D* f0 = faces[0];
					Face3D* f1 = faces[1];
					Face3D* f2 = faces[2];
					Face3D* f3 = faces[3];
#endif

					cout << "BAD TET KEY: " << (int)key << endl;
					exit(-1);
				}

				int new_index[11];
				if(parity_flip[tet->tet_index])
				{
					for(int i=0; i < 11; i++)
						new_index[i] = vertexTableEven[key][i+4];
				}
				else
				{
					for(int i=0; i < 11; i++)
						new_index[i] = vertexTableOdd[key][i+4];
				}

				Vertex3D *v[15];
				for(int i=0; i < 4; i++)
					v[i] = verts[i];

				v[4] = edges[0]->cut;
				v[5] = edges[1]->cut;
				v[6] = edges[2]->cut;
				v[7] = edges[3]->cut;
				v[8] = edges[4]->cut;
				v[9] = edges[5]->cut;
				v[10] = faces[0]->triple;
				v[11] = faces[1]->triple;
				v[12] = faces[2]->triple;
				v[13] = faces[3]->triple;
				v[14] = tet->quad;

				// before applying changes, we can check if we agree existing decision
				// if no quad, some cuts are also missing
				for(int e=0; e < 6; e++)
				{
					// debug check
					if(edges[e]->cut != NULL)
					{
						if(edges[e]->cut != v[new_index[e]])
							cout << "Warning: Cut Generalization Disagrees with neighbor tet!" << endl;
					}
					edges[e]->cut = v[new_index[e]];
				}

				// if no quad, some triples are also missing
				for(int f=0; f < 4; f++)
				{
					// debug check
					if(faces[f]->triple != NULL)
					{
						if(faces[f]->triple != v[new_index[f+6]])
							cout << "Warning: Face Generalization Disagrees with neighbor tet!" << endl;
					}
					faces[f]->triple = v[new_index[f+6]];
				}

				// set quad, which we know is missing
				tet->quad = v[new_index[_ABCD-4]];

				// verify that this tet is now generalized
				Vertex3D *v2[15];
				memset(v2, 0, 15*sizeof(Vertex3D*));
				lattice->getRightHandedVertexList(tet, v2);
				for(int i=0; i < 15; i++)
				{
					if(v2[i] == NULL)
					{
						cout << "Problem! Tet failed to generalize" << endl;
						exit(-1);
					}
				}

			}
			// end generalization
			else
			{
				// set key for non-degenerate case
				Edge3D *test_edges[EDGES_PER_TET];
				lattice->getEdgesAroundTet(tet, test_edges);
				tet->key = lattice->keyFromAdjacentEdges(test_edges);
			}
		}
	}


	//--------------------------------------------------
	//      Generalize All Tets In Buffer Cells
	//--------------------------------------------------
	for(unsigned int c = 0; c < lattice->buffer_cells.size(); c++)
	{

		OTCell *cell = lattice->buffer_cells[c];

		for(int t=0; t < TETS_PER_CELL; t++)
		{
			Tet3D *tet = cell->tets[t];
			if(!tet)
				continue;

			// if no quad, start generalization
			if(!tet->quad)
			{
				// look up generalization
				Vertex3D *verts[4];
				Edge3D *edges[6];
				Face3D *faces[4];

				lattice->getAdjacencyLists(tet, verts, edges, faces);
				unsigned char key = tet->key = lattice->keyFromAdjacentEdges(edges);

				if(!lattice->isKeyValid(key)){
					cout << "BAD TET KEY: " << key << endl;
					exit(-1);
				}

				int new_index[11];
				if(parity_flip[tet->tet_index])
				{
					for(int i=0; i < 11; i++)
						new_index[i] = vertexTableEven[key][i+4];
				}
				else
				{
					for(int i=0; i < 11; i++)
						new_index[i] = vertexTableOdd[key][i+4];
				}

				Vertex3D *v[15];
				for(int i=0; i < 4; i++)
					v[i] = verts[i];

				v[4] = edges[0]->cut;
				v[5] = edges[1]->cut;
				v[6] = edges[2]->cut;
				v[7] = edges[3]->cut;
				v[8] = edges[4]->cut;
				v[9] = edges[5]->cut;
				v[10] = faces[0]->triple;
				v[11] = faces[1]->triple;
				v[12] = faces[2]->triple;
				v[13] = faces[3]->triple;
				v[14] = tet->quad;

				// before applying changes, we can check if we agree existing decision
				// if no quad, some cuts are also missing
				for(int e=0; e < 6; e++){
					// debug check
					if(edges[e]->cut != NULL)
					{
						if(edges[e]->cut != v[new_index[e]])
							cout << "Warning: Cut Generalization Disagrees with neighbor tet!" << endl;
					}
					edges[e]->cut = v[new_index[e]];
				}

				// if no quad, some triples are also missing
				for(int f=0; f < 4; f++){
					// debug check
					if(faces[f]->triple != NULL)
					{
						if(faces[f]->triple != v[new_index[f+6]])
							cout << "Warning: Face Generalization Disagrees with neighbor tet!" << endl;
					}
					faces[f]->triple = v[new_index[f+6]];
				}

				// set quad, which we know is missing
				tet->quad = v[new_index[_ABCD-4]];

				// verify that this tet is now generalized
				Vertex3D *v2[15];
				memset(v2, 0, 15*sizeof(Vertex3D*));
				lattice->getRightHandedVertexList(tet, v2);
				for(int i=0; i < 15; i++)
				{
					if(v2[i] == NULL)
					{
						cout << "Unhandled Problem! Tet failed to generalize" << endl;
						exit(-1);
					}
				}

			}
			// end generalization
			else
			{
				// set key for non-degenerate case
				Edge3D *test_edges[EDGES_PER_TET];
				lattice->getEdgesAroundTet(tet, test_edges);
				tet->key = lattice->keyFromAdjacentEdges(test_edges);
			}
		}
	}

	lattice->setGeneralized(true);
}

void Cleaver::BCCLattice3DMesherLite::fill_all_stencils()
{
	//--------------------------------------
	//      fill cut stencils
	//--------------------------------------
	for(unsigned int c = 0; c < lattice->cut_cells.size(); c++)
	{
		OTCell *cell = lattice->cut_cells[c];

		for(int t=0; t < TETS_PER_CELL; t++)
		{
			Tet3D *tet = cell->tets[t];
			if(!tet->stenciled)
				fill_stencil(tet);
		}
	}

	//--------------------------------------
	//      fill buffer stencils
	//--------------------------------------
	for(unsigned int c = 0; c < lattice->buffer_cells.size(); c++)
	{
		OTCell *cell = lattice->buffer_cells[c];

		for(int t=0; t < TETS_PER_CELL; t++)
		{
			Tet3D *tet = cell->tets[t];
			if(tet && !tet->stenciled)
				fill_stencil(tet);
		}
	}

	lattice->setStenciled(true);
}

void Cleaver::BCCLattice3DMesherLite::compute_cut( Edge3D *edge )
{
	Vertex3D *v1 = edge->v1;
	Vertex3D *v2 = edge->v2;

#ifdef _DEBUG

	vec3 v1pos = v1->pos();
	vec3 v2pos = v2->pos();
	if (  length(v1pos - vec3(24.5,8.5,6.5)) < 1e-7 && length(v2pos - vec3(24.0,8.0,6.0)) < 1e-7 )
	{
		int aaa = 0;
	}

#endif
	

	edge->evaluated = true;
	if(!isTransition(v1, v2))
		return;

	int a_mat = v1->label;
	int b_mat = v2->label;

	double t = 0.5;
	Vertex3D *cut = new Vertex3D(lattice->materials());    
	cut->pos() = v1->pos()*(1-t) + v2->pos()*t;

	cut->closestGeometry = v1;

	//cut->vals[a_mat] = mat_val;
	//cut->vals[b_mat] = mat_val;
	cut->label = a_mat;         // doesn't really matter which
	//cut->lbls[v1->label] = true;
	//cut->lbls[v2->label] = true;
	cut->violating = false;

	//lattice->cuts.push_back(cut);
	edge->cut = cut;
	//cut->m_edge = edge;
	cut->order() = 1;
}

bool Cleaver::BCCLattice3DMesherLite::isTransition( const Vertex3D* v1, const Vertex3D* v2 )
{
	int v1_label = v1->label;
	int v2_label = v2->label;

	return (v1_label != v2_label);
}

void Cleaver::BCCLattice3DMesherLite::fixTriangleOrdering( Edge3D *edges[], Vertex3D *verts[] )
{
	Edge3D *tmp_edge;

	// fix e[0] to match v[0]
	for(int e=0; e < 3; e++)
	{
		if(edges[e]->v1 != verts[0] &&
			edges[e]->v2 != verts[0])
		{
			// swap
			tmp_edge = edges[0];
			edges[0] = edges[e];
			edges[e] = tmp_edge;
		}
	}

	// fix e[1] to match v[1]
	for(int e=1; e < 3; e++)
	{
		if(edges[e]->v1 != verts[1] &&
			edges[e]->v2 != verts[1])
		{
			// swap
			tmp_edge = edges[1];
			edges[1] = edges[e];
			edges[e] = tmp_edge;
		}
	}
}

void Cleaver::BCCLattice3DMesherLite::compute_triple( Face3D *face )
{
	Vertex3D *verts[3];
	Edge3D *edges[3];

	lattice->getAdjacencyLists(face, verts, edges);

	face->evaluated = true;
	if(!edges[0]->cut || !edges[1]->cut || !edges[2]->cut) // any one cut is NULL, then no triple. 
		return;

	fixTriangleOrdering(edges, verts);

	//-------------------------------------------------------
	// get coordinates and create Vertex
	//-------------------------------------------------------
	Vertex3D* v1 = verts[0];
	Vertex3D* v2 = verts[1];
	Vertex3D* v3 = verts[2];

	vec3 result = vec3::zero;

	float a1, b1, c1, d1;
	float a2, b2, c2, d2;
	float a3, b3, c3, d3;

	// get materials
	int m1 = v1->label;
	int m2 = v2->label;
	int m3 = v3->label;

	int axis = LATTICE_FACE_AXIS[face->face_index];

	result = (1.0/3.0)*(v1->pos() + v2->pos() + v3->pos());

	Vertex3D *triple = new Vertex3D(lattice->materials());
	triple->pos() = result;
	//triple->lbls[verts[0]->label] = true;
	//triple->lbls[verts[1]->label] = true;
	//triple->lbls[verts[2]->label] = true;

	triple->order() = TRIP;
	triple->violating = false;
	triple->closestGeometry = NULL;    
	face->triple = triple;

	// check if its violating
	return;
}

void Cleaver::BCCLattice3DMesherLite::compute_quadruple( Tet3D *tet )
{
	Vertex3D *verts[4];
	Edge3D *edges[6];
	Face3D *faces[4];

	lattice->getAdjacencyLists(tet, verts, edges, faces);

	fixTetrahedronOrdering(faces, edges, verts);

	tet->evaluated = true;
	if(!(edges[0]->cut && edges[0]->cut->order() == 1 &&
		edges[1]->cut && edges[1]->cut->order() == 1 &&
		edges[2]->cut && edges[2]->cut->order() == 1 &&
		edges[3]->cut && edges[3]->cut->order() == 1 &&
		edges[4]->cut && edges[4]->cut->order() == 1 &&
		edges[5]->cut && edges[5]->cut->order() == 1))
		return;

	Vertex3D* v1 = verts[0];
	Vertex3D* v2 = verts[1];
	Vertex3D* v3 = verts[2];
	Vertex3D* v4 = verts[3];

	vec3 result = vec3::zero;
	result = (1.0/4.0)*(v1->pos() + v2->pos() + v3->pos() + v4->pos());

	Vertex3D *quad = new Vertex3D(lattice->materials());
	quad->pos() = result;
	//quad->lbls[verts[0]->label] = true;
	//quad->lbls[verts[1]->label] = true;
	//quad->lbls[verts[2]->label] = true;
	//quad->lbls[verts[3]->label] = true;

	quad->order() = QUAD;
	quad->violating = false;
	quad->closestGeometry = NULL;
	tet->quad = quad;
	return;
}

void Cleaver::BCCLattice3DMesherLite::fixTetrahedronOrdering( Face3D *faces[], Edge3D *edges[], Vertex3D *verts[] )
{
	for(int j=0; j < 3; j++)
	{
		// fix f[j] to match v[j]
		for(int f=j; f < 4; f++)
		{
			Vertex3D *tri_verts[3];
			lattice->getVertsAroundFace(faces[f], tri_verts);

			// if no vertices shared, must be opposite
			if(tri_verts[0] != verts[j] &&
				tri_verts[1] != verts[j] &&
				tri_verts[2] != verts[j])
			{
				// swap
				Face3D *tmp_face = faces[j];
				faces[j] = faces[f];
				faces[f] = tmp_face;
			}
		}
	}
}

void Cleaver::BCCLattice3DMesherLite::fill_stencil( Tet3D *tet )
{
	Vertex3D *verts[15];
	lattice->getRightHandedVertexList(tet, verts);
	unsigned char key = 63; // tet->key;

	if(parity_flip[tet->tet_index])
	{//odd
		for(int t=0; t < 24; t++) //every stencil contains at most 24 tets
		{
			//--------------------------------------
			//  Proceed Only If Should Output Tets
			//--------------------------------------
			if(stencilTableOdd[key][t][0] == _O)
				break;

			//-----------------------
			//     Get Vertices
			//-----------------------
			Vertex3D *v1 = verts[stencilTableOdd[key][t][0]]->root();  // grabbing root ensures uniqueness
			Vertex3D *v2 = verts[stencilTableOdd[key][t][1]]->root();
			Vertex3D *v3 = verts[stencilTableOdd[key][t][2]]->root();
			Vertex3D *v4 = verts[stencilTableOdd[key][t][3]]->root();
			Vertex3D *vM = verts[materialTableOdd[key][t]]->root();

			//-----------------------
			//  Check If Degenerate
			//-----------------------
			if(v1 == v2 || v1 == v3 || v1 == v4 || v2 == v3 || v2 == v4 || v3 == v4)
				continue;

			//----------------------------
			//  Create Tet + Add to List
			//----------------------------
			Tet *st = lattice->tree->createTet(v1, v2, v3, v4, (int)vM->label);
			st->key = tet->key;
		}
	}
	else 
	{//even
		for(int t=0; t < 24; t++)
		{
			//--------------------------------------
			//  Procede Only If Should Output Tets
			//--------------------------------------
			if(stencilTableEven[key][t][0] == _O)
				break;

			//-----------------------
			//     Get Vertices
			//-----------------------
			Vertex3D *v1 = verts[stencilTableEven[key][t][0]]->root();
			Vertex3D *v2 = verts[stencilTableEven[key][t][1]]->root();
			Vertex3D *v3 = verts[stencilTableEven[key][t][2]]->root();
			Vertex3D *v4 = verts[stencilTableEven[key][t][3]]->root();
			Vertex3D *vM = verts[materialTableEven[key][t]]->root();

			//-----------------------
			//  Check If Degenerate
			//-----------------------
			if(v1 == v2 || v1 == v3 || v1 == v4 || v2 == v3 || v2 == v4 || v3 == v4)
				continue;

			//----------------------------
			//  Create Tet + Add to List
			//----------------------------
			Tet *st = lattice->tree->createTet(v1, v2, v3, v4, vM->label);
			st->key = tet->key;
		}

	}


	tet->stenciled = true;
}
