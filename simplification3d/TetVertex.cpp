#include "stdafx.h"
#include "TetVertex.h"
#include "TetrahedralMesh.h"

using namespace std;
using namespace cv;

int cwg::TetVertex::ms_nextid = 0;
extern SIMPLIFICATION3D_API float glb_theta;

void cwg::TetVertex::ComputeInitialQ()
{
	m_Q.set2zero();
	//m_Q.diagval = 1.0;
	m_Q(0,0) = 1.0;
	m_Q(1,1) = 1.0;
	m_Q(2,2) = 1.0;
	m_Q(0,3) = -m_xyz.x;
	m_Q(1,3) = -m_xyz.y;
	m_Q(2,3) = -m_xyz.z;
	m_Q(3,3) = m_xyz.dot(m_xyz);

}

void cwg::TetVertex::ComputeInitialQ(TetrahedralMesh* tetmesh)
{
	m_Q.set2zero();
	//m_Q.diagval = 1.0;
	m_Q(0,0) = 1.0;
	m_Q(1,1) = 1.0;
	m_Q(2,2) = 1.0;
	m_Q(0,3) = -m_xyz.x;
	m_Q(1,3) = -m_xyz.y;
	m_Q(2,3) = -m_xyz.z;
	m_Q(3,3) = m_xyz.dot(m_xyz);

	SymMat4 tmpQ;
	tmpQ.set2zero();

	if (this->NLabels() != 1)
	{		
		avg_weight = 100;
		for (int i=0; i<OrigBoundaryTris.size(); i++)
		{
			TetTriangle* tri = tetmesh->InternalMesh->GetTri( OrigBoundaryTris[i] );

			SymMat4 Qi = compute_tri_Q(&tetmesh->m_storage_verts[ tri->Vert(0) ], 
				&tetmesh->m_storage_verts[ tri->Vert(1) ], 
				&tetmesh->m_storage_verts[ tri->Vert(2) ]);
			tmpQ += Qi;


		}
		m_Q = (m_Q + tmpQ);
	}
	
	m_Q *= avg_weight;
}


cwg::SymMat4 cwg::compute_tri_Q( TetVertex* v0, TetVertex* v1, TetVertex* v2 )
{
	Cleaver::vec3 v0vec3 = v0->pos();
	Cleaver::vec3 v1vec3 = v1->pos();
	Cleaver::vec3 v2vec3 = v2->pos();
	Cleaver::vec3 v0v1 = v1vec3 - v0vec3;
	Cleaver::vec3 v0v2 = v2vec3 - v0vec3;
	Cleaver::vec3 n = v0v1.cross(v0v2);

	double d = (-n.x*v0->pos().x -n.y*v0->pos().y -n.z*v0->pos().z)/length(n);

	n = Cleaver::normalize(n);
	double g = n.dot(v0vec3);
	Cleaver::vec4 nm(n.x,n.y, n.z, d);
	SymMat4 Q = nm.vvt();
	return Q;
}



void cwg::TetVertex::invalidate()
{
	ID = -1;
	clear();
}

void cwg::TetVertex::update_labels( TetrahedralMesh* tetmesh )
{
	m_labels.clear();
	for (vector<int>::iterator it = Tets.begin(); it != Tets.end(); it++)
	{
		int lb = tetmesh->m_storage_tets[ *it ].get_label();
		assign_label(lb);
		
	}
}

void cwg::TetVertex::read( std::ifstream& ifile )
{
	ifile.read( (char*)&ID, sizeof(int) );
	double vals[5];
	ifile.read( (char*)vals, sizeof(double)*5 );
	m_xyz.x = vals[0];
	m_xyz.y = vals[1];
	m_xyz.z = vals[2];
	m_funcval = vals[3];
	setAvgWeight(vals[4]);
	ifile.read( (char*)&m_border_code, sizeof(unsigned char) );
}

void cwg::TetVertex::write( std::ofstream& ofile ) const
{
	ofile.write( (char*)&ID, sizeof(int) );
	double vals[5];
	vals[0] = m_xyz.x;
	vals[1] = m_xyz.y;
	vals[2] = m_xyz.z;
	vals[3] = m_funcval;
	vals[4] = avg_weight;
	ofile.write( (char*)vals, sizeof(double)*5 );
	ofile.write( (char*)&m_border_code, sizeof(unsigned char) );
}

void cwg::TetVertex::copy_info_from_cleaver_vert( Cleaver::Vertex3D* v )
{
	m_xyz = v->pos();
	m_funcval = m_xyz.dot(m_xyz);
}

void cwg::TetVertex::clear()
{
	Tets.clear();
	//m_Q.release();
	m_labels.clear();
}

void cwg::TetVertex::correlate_tet( int tet )
{
	//Tets.insert(tet);
	std::vector<int>::iterator it = std::lower_bound(Tets.begin(), Tets.end(), tet);
	if (it != Tets.end())
	{
		if (*it != tet)
			Tets.insert(it, tet);
	}
	else
		Tets.push_back(tet);
}

void cwg::TetVertex::decorrelate_tet( int tet )
{
	//Tets.erase(tet);
	std::vector<int>::iterator it = std::lower_bound(Tets.begin(), Tets.end(), tet);
	if (it != Tets.end())
	{
		if (*it != tet)
 			throw exception("void cwg::TetVertex::decorrelate_tet( int tet ): tet not exists in the vert! (if-if) === Slide window is too small or simplify too many");
		else
			Tets.erase(it);
	}
	else
		throw exception("void cwg::TetVertex::decorrelate_tet( int tet ): tet not exists in the vert! (else) === Slide window is too small or simplify too many");
	
}

void cwg::TetVertex::correlate_boundary_tri( int tri )
{
	std::vector<int>::iterator it = std::lower_bound(OrigBoundaryTris.begin(), OrigBoundaryTris.end(), tri);
	if (it != OrigBoundaryTris.end())
	{
		if (*it != tri)
			OrigBoundaryTris.insert(it, tri);
	}
	else
		OrigBoundaryTris.push_back(tri);
}

void cwg::TetVertex::assign_label(int lb )
{
	std::vector<int>::iterator it = std::lower_bound(m_labels.begin(), m_labels.end(), lb);
	if (it != m_labels.end())
	{
		if (*it != lb)
			m_labels.insert(it, lb);
	}
	else
		m_labels.push_back(lb);
}

std::istream& cwg::operator>>( std::istream& is, cwg::TetVertex& obj )
{
	double x, y, z;
	is >> obj.ID >> x >> y >> z >> obj.m_funcval >> obj.avg_weight >> obj.m_border_code; 
	obj.m_xyz = Cleaver::vec3(x, y, z);
	return is;
}

std::ostream& cwg::operator<<( std::ostream& os, const TetVertex& obj )
{
	os << obj.ID << " " << obj.m_xyz << " " << obj.m_funcval << " " << obj.avg_weight << " " << (int)obj.m_border_code;
	return os;
}