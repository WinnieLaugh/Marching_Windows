#ifndef _CWG_BASIC_GEOMETRY_H_
#define _CWG_BASIC_GEOMETRY_H_

#include <cstring>
#include <vector>
#include "vec3.h"
#include "Vertex.h"
#include "Octree.h"

namespace Cleaver
{
	class Edge3D : public Geometry
	{

	public:
		Edge3D():v1(0),v2(0),cut(0),isLong(false),evaluated(false){}
		Edge3D(OTCell *cell, int index):cell(cell), v1(0), v2(0), cut(0), edge_index(index), isLong(false), evaluated(false){}
		Edge3D(bool isLong, OTCell *cell, int index):cell(cell), v1(0), v2(0), cut(0), edge_index(index), isLong(isLong), evaluated(false){}
		virtual ~Edge3D();

		OTCell *cell;    
		Vertex3D *v1, *v2;      // two ordered vertices adjacent to edge
		Vertex3D *cut;

		unsigned char edge_index;
		bool isLong:1;
		bool evaluated:1;

		float length()
		{
			return (float)L2(v1->pos() - v2->pos());
		}

		bool contains(Vertex3D *vertex)
		{
			if(this->v1 == vertex || this->v2 == vertex)
				return true;
			else
				return false;
		}

		bool containsBoth(Vertex3D *v1, Vertex3D *v2)
		{
			if((this->v1 == v1 || this->v2 == v1) &&
				(this->v1 == v2 || this->v2 == v2))
				return true;
			else
				return false;
		}


	};

	class Face3D : public Geometry
	{
	public:
		Face3D(Vertex3D *t) : triple(t),evaluated(true){ }
		Face3D():triple(NULL),evaluated(false){ }
		Face3D(OTCell *cell, int index):cell(cell), triple(NULL), face_index(index), evaluated(false){ }
		~Face3D();

		OTCell *cell;
		Vertex3D *triple;
		unsigned char face_index;
		bool evaluated:1;
	};

	class Tet3D : public Geometry
	{
	public:
		Tet3D():quad(NULL),tm_index(0),key(0),evaluated(false),stenciled(false){}
		Tet3D(int index):quad(NULL),tm_index(0),tet_index(index),key(0),evaluated(false),stenciled(false){}
		Tet3D(OTCell *cell,int index):cell(cell),quad(NULL), tm_index(0),tet_index(index),key(0),evaluated(false),stenciled(false){}
		~Tet3D();

		OTCell *cell;    
		Vertex3D *quad;

		int *tm_index;
		unsigned char tet_index;
		unsigned char key;

		bool evaluated:1;
		bool stenciled:1;
	};


#define DATA3D(a,b,c,m,mat) data[m*(a + b*w + c*w*h) + mat]
#define LBL3D(a,b,c)    labels[(a + b*w + c*w*h)]
#define LBL3DOCTREE(a,b,c)    labels[(a + b*(w+1) + c*(w+1)*(h+1))]
#define DUAL3D(a,b,c)   dual[a + b*(w-1) + c*(w-1)*(h-1)]

#define MAIN_HORI(a,b,c) main_hori[a + b*(w-1) + c*(w-1)*h]
#define MAIN_VERT(a,b,c) main_vert[a + b*w + c*w*(h-1)]
#define MAIN_DEEP(a,b,c) main_deep[a + b*w + c*(w-1)*(h-1)]

#define DUAL_HORI(a,b,c) dual_hori[a + b*(w-1) + c*(w-1)*h]
#define DUAL_VERT(a,b,c) dual_vert[a + b*w + c*w*(h-1)]
#define DUAL_DEEP(a,b,c) dual_deep[a + b*w + c*(w-1)*(h-1)]


#define DIAG_ULF(a,b,c) diagULF[a + b*(w-1) + c*(w-1)*(h-1)]
#define DIAG_ULB(a,b,c) diagULB[a + b*(w-1) + c*(w-1)*(h-1)]
#define DIAG_URF(a,b,c) diagURF[a + b*(w-1) + c*(w-1)*(h-1)]
#define DIAG_URB(a,b,c) diagURB[a + b*(w-1) + c*(w-1)*(h-1)]
#define DIAG_LLF(a,b,c) diagLLF[a + b*(w-1) + c*(w-1)*(h-1)]
#define DIAG_LLB(a,b,c) diagLLB[a + b*(w-1) + c*(w-1)*(h-1)]
#define DIAG_LRF(a,b,c) diagLRF[a + b*(w-1) + c*(w-1)*(h-1)]
#define DIAG_LRB(a,b,c) diagLRB[a + b*(w-1) + c*(w-1)*(h-1)]
}

#endif