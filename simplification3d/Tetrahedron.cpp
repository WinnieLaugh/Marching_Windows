#include "stdafx.h"
#include "Tetrahedron.h"
#include "TetrahedralMesh.h"
#include "Utils/quadratic_integral_tet.h"

using namespace std;
using namespace cv;

int cwg::Tetrahedron::ms_nextid = 0;

void cwg::Tetrahedron::correlate_verts(TetrahedralMesh* tetmesh)
{
	for (int i=0; i<4; i++)
	{
		tetmesh->m_storage_verts[ this->Verts[i] ].assign_label(m_label);
		tetmesh->m_storage_verts[ this->Verts[i] ].correlate_tet(ID);
	}
}
void cwg::Tetrahedron::correlate_verts_sliding_result(TetrahedralMesh* tetmesh)
{
	for (int i=0; i<4; i++)
	{
		tetmesh->ResultV[ this->Verts[i] ].assign_label(m_label);
		tetmesh->ResultV[ this->Verts[i] ].correlate_tet(ID);
	}
}

int cwg::Tetrahedron::find_vert( int vid ) const
{
	for (int i=0; i<4; i++)
	{
		if (Verts[i] == vid)
			return i;
	}
	return -1;
}


void cwg::Tetrahedron::invalidate(TetrahedralMesh* tetmesh)
{
	for(int index=0; index<4; index++)
		Verts[index] = -1;

	ID = -1;
	m_label = -1;

}

bool cwg::Tetrahedron::inside_tet_check(float x, float y, float z, std::vector<TetVertex> *vertex_vector){
	Cleaver::vec3 p;
	p.x = x;
	p.y = y;
	p.z = z;

	Cleaver::vec3 a = (*vertex_vector)[ this->Verts[0] ].pos();
	Cleaver::vec3 b = (*vertex_vector)[ this->Verts[1] ].pos();
	Cleaver::vec3 c = (*vertex_vector)[ this->Verts[2] ].pos();
	Cleaver::vec3 d = (*vertex_vector)[ this->Verts[3] ].pos();

	float v0 = cwg::tet_vol(a, b, c, d);
	float v1 = cwg::tet_vol(p, b, c, d);
	float v2 = cwg::tet_vol(a, p, c, d);
	float v3 = cwg::tet_vol(a, b, p, d);
	float v4 = cwg::tet_vol(a, b, c, p);

	if (((v0 > 0) && (v1 > 0) && (v2 > 0) && (v3 > 0) && (v4 > 0))||
		((v0 < 0) && (v1 < 0) && (v2 < 0) && (v3 < 0) && (v4 < 0)))
	{
		return true;
	}else{
		return false;
	}

}

double cwg::Tetrahedron::variance( const TetrahedralMesh* tetmesh , double avg) const
{
	Cleaver::vec3 a = tetmesh->m_storage_verts[ this->Verts[0] ].pos();
	Cleaver::vec3 b = tetmesh->m_storage_verts[ this->Verts[1] ].pos();
	Cleaver::vec3 c = tetmesh->m_storage_verts[ this->Verts[2] ].pos();
	Cleaver::vec3 d = tetmesh->m_storage_verts[ this->Verts[3] ].pos();

	double e[6];
	e[0] = length(a-b);
	e[1] = length(a-c);
	e[2] = length(a-d);
	e[3] = length(b-c);
	e[4] = length(b-d);
	e[5] = length(c-d);

	double sum = 0;
	for (int i=0; i<6; ++i)
	{
		sum += e[i];
	}
	avg = sum/6;
	double v=0;
	for (int i=0; i<6; ++i)
	{
		v += (e[i]-avg)*(e[i]-avg);
	}

	return v/6;
	
}

double cwg::Tetrahedron::avg_edge_length_sliding_result( const TetrahedralMesh* tetmesh) const
{
	Cleaver::vec3 a = tetmesh->ResultV[ this->Verts[0] ].pos();
	Cleaver::vec3 b = tetmesh->ResultV[ this->Verts[1] ].pos();
	Cleaver::vec3 c = tetmesh->ResultV[ this->Verts[2] ].pos();
	Cleaver::vec3 d = tetmesh->ResultV[ this->Verts[3] ].pos();

	double e[6];
	e[0] = length(a-b);
	e[1] = length(a-c);
	e[2] = length(a-d);
	e[3] = length(b-c);
	e[4] = length(b-d);
	e[5] = length(c-d);

	double sum = 0;
	for (int i=0; i<6; ++i)
	{
		sum += e[i];
	}
	
	return sum/6;
	
}


double cwg::Tetrahedron::compute_volume_sliding_result( const TetrahedralMesh* tetmesh ) const
{

	Cleaver::vec3 a = tetmesh->ResultV[ this->Verts[0] ].pos();
	Cleaver::vec3 b = tetmesh->ResultV[ this->Verts[1] ].pos();
	Cleaver::vec3 c = tetmesh->ResultV[ this->Verts[2] ].pos();
	Cleaver::vec3 d = tetmesh->ResultV[ this->Verts[3] ].pos();
	
	return cwg::tet_vol(a, b, c, d);
}
double cwg::Tetrahedron::compute_volume( const TetrahedralMesh* tetmesh ) const
{
	if(this->Verts[0] == -1 && this->Verts[1] == -1 && this->Verts[2] == -1 && this->Verts[3] == -1){
		std::cout << "where the error begins" << std::endl;
		std::cout << "ID " << this->ID << std::endl;

		std::cout << "Verts: " << this->Verts[0] << " " << this->Verts[1] << " " << this->Verts[2] << " " << this->Verts[3] << std::endl;
		std::cout << "=============================test==========================" << std::endl;
		std::cout << "m_storage_verts Vert 0 x: " << tetmesh->m_storage_verts[ this->Verts[0] ].pos().x << std::endl;
	}

	Cleaver::vec3 a = tetmesh->m_storage_verts[ this->Verts[0] ].pos();
	Cleaver::vec3 b = tetmesh->m_storage_verts[ this->Verts[1] ].pos();
	Cleaver::vec3 c = tetmesh->m_storage_verts[ this->Verts[2] ].pos();
	Cleaver::vec3 d = tetmesh->m_storage_verts[ this->Verts[3] ].pos();

	return cwg::tet_vol(a, b, c, d);
}

void cwg::Tetrahedron::read( std::ifstream& ifile )
{
	ifile.read( (char*)&ID, sizeof(int) );
	ifile.read( (char*)Verts, sizeof(int)*4 );
	ifile.read( (char*)&m_label, sizeof(int) );
}

void cwg::Tetrahedron::write( std::ofstream& ofile ) const
{
	ofile.write( (char*)&ID, sizeof(int) );
	ofile.write( (char*)Verts, sizeof(int)*4 );
	ofile.write( (char*)&m_label, sizeof(int) );
}

std::vector<double> cwg::Tetrahedron::compute_dihedral_angle( const TetrahedralMesh* tetmesh ) const
{
	vector<TetVertex*> verts(4);
	vector<Cleaver::vec3> verts_pos(4);
	for (int i=0; i<4; i++) 
	{
		verts_pos[i] = tetmesh->m_storage_verts[ this->Verts[i] ].pos();
	}
	return cwg::tet_dihedral_angles(verts_pos[0], verts_pos[1], verts_pos[2], verts_pos[3]);
}

double cwg::Tetrahedron::minimal_dihedral_angle( const TetrahedralMesh* tetmesh ) const
{
	vector<double> angles = compute_dihedral_angle(tetmesh);
	double mda = numeric_limits<double>::max();
	for (int i=0; i<angles.size(); i++)
		mda = std::min(mda, angles[i]);
	return mda;
}

void cwg::Tetrahedron::compute_tri_ratio(TetrahedralMesh* tetmesh, std::vector<double>& surf_ratio) const
{
	double R,r;
	
	for (int i=0; i<tetmesh->Tris.size(); i++)
	{
		int t0 = tetmesh->Tris[i]->Tets[0];
		int t1 = tetmesh->Tris[i]->Tets[1];

		bool is_a_surface = false;
		
		if (t0 < 0 && t1 < 0)
			throw exception("bool cwg::TetrahedralMesh::save_txt_surf_ply( const std::string& filename ) const:\
							t0 < 0 && t1 < 0");
		else if (t0 >= 0 && t1 < 0)
		{
			is_a_surface = true;
			int lb0 = tetmesh->m_storage_tets[t0].get_label();
			//if(lb0 > 0)
			{
				R = tri_circ_r(tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(0)].pos(), tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(1)].pos()
					,tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(2)].pos());
				r = tri_in_r(tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(0)].pos(), tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(1)].pos()
					,tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(2)].pos());
				surf_ratio.push_back(2*r/R);
			}
			
		}
		else if (t1 >= 0 && t0 < 0)
		{
			is_a_surface = true;
			int lb1 = tetmesh->m_storage_tets[t1].get_label();
			//if(lb1 > 0)
			{
				R = tri_circ_r(tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(0)].pos(), tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(1)].pos()
					,tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(2)].pos());
				r = tri_in_r(tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(0)].pos(), tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(1)].pos()
					,tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(2)].pos());
				surf_ratio.push_back(2*r/R);
			}
		}
		else
		{
			int lb0 = tetmesh->m_storage_tets[t0].get_label();
			int lb1 = tetmesh->m_storage_tets[t1].get_label();
			if (lb0 != lb1)
			{
				is_a_surface = true;
				
				R = tri_circ_r(tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(0)].pos(), tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(1)].pos()
					,tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(2)].pos());
				r = tri_in_r(tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(0)].pos(), tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(1)].pos()
					,tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(2)].pos());
				surf_ratio.push_back(2*r/R);
				
			}
		}

		
	}
	
	
}

void cwg::Tetrahedron::compute_tri_ratio_cgal(TetrahedralMesh* tetmesh, std::vector<double>& surf_ratio) const
{
	double R,r;

	for (int i=0; i<tetmesh->Tris.size(); i++)
	{
		
		R = tri_circ_r(tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(0)].pos(), tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(1)].pos()
			,tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(2)].pos());
		r = tri_in_r(tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(0)].pos(), tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(1)].pos()
			,tetmesh->m_storage_verts[tetmesh->Tris[i]->Vert(2)].pos());
		surf_ratio.push_back(2*r/R);
	}


}
double cwg::tri_circ_r( const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
					   const Cleaver::vec3& v2pos )
{
	double area = tri_area(v0pos, v1pos, v2pos);
	double a = length( v2pos - v1pos);
	double b = length(v1pos - v0pos);
	double c = length(v2pos - v0pos);

	double	r = a*b*c / (4*area);
	return r;
}

double cwg::tri_in_r( const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
					 const Cleaver::vec3& v2pos )
{
	double area = tri_area(v0pos, v1pos, v2pos);
	double a = length( v2pos - v1pos);
	double b = length(v1pos - v0pos);
	double c = length(v2pos - v0pos);
	

	double r = 2 * area / (a + b + c);
	return r;
}
double cwg::tri_area( const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
					 const Cleaver::vec3& v2pos )
{
	return length(cross(v1pos - v0pos, v2pos - v0pos)) * 0.5;
}

Cleaver::vec3& cwg::tri_circ_c( const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
							   const Cleaver::vec3& v2pos )
{
	Cleaver::vec3 e1 = v1pos - v0pos;
	Cleaver::vec3 e2 = v2pos - v0pos;

	Cleaver::vec3 tNorm = cross(e1,e2);

	Cleaver::vec3 e1Norm = cross(e1,tNorm);
	Cleaver::vec3 e2Norm = cross(e2,tNorm);

	Cleaver::vec3 midv1 = (v1pos + v0pos) / 2;
	Cleaver::vec3 midv2 = (v2pos + v0pos) / 2;

	//L1: x=x1+nx1*t, y=y1+ny1*t,z=z1+nz1*t 
	//L2: x=x2+nx2*s, y=y2+ny2*s,z=z2+nz2*s
	//find intersection
	double t = ((midv1.y - midv2.y)*e2Norm.x - e2Norm.y * (midv1.x - midv2.x))/(e2Norm.y * e1Norm.x - e2Norm.x * e1Norm.y);
	//double s = (midv1.x - midv2.x + e1Norm.x * t)/ e2Norm.x;
	
	//Cleaver::vec3 circ_c = midv1 + e1Norm * t;

	return midv1 + e1Norm * t;
}

void cwg::Tetrahedron::compute_sphere_ratio( TetrahedralMesh* tetmesh ) const
{
	double R,r;
	for (int i=0; i<tetmesh->m_storage_tets.size(); i++)
	{
		R = tet_circ_r(tetmesh->m_storage_verts[tetmesh->m_storage_tets[i].Verts[0]].pos(), tetmesh->m_storage_verts[tetmesh->m_storage_tets[i].Verts[1]].pos(),
			tetmesh->m_storage_verts[tetmesh->m_storage_tets[i].Verts[2]].pos(), tetmesh->m_storage_verts[tetmesh->m_storage_tets[i].Verts[3]].pos());
		r = tet_in_r(tetmesh->m_storage_verts[tetmesh->m_storage_tets[i].Verts[0]].pos(), tetmesh->m_storage_verts[tetmesh->m_storage_tets[i].Verts[1]].pos(),
			tetmesh->m_storage_verts[tetmesh->m_storage_tets[i].Verts[2]].pos(), tetmesh->m_storage_verts[tetmesh->m_storage_tets[i].Verts[3]].pos());

		tetmesh->m_storage_tets[i].set_ratio(3*r/R);
	}
	
}
double cwg::Tetrahedron::compute_tet_ratio( TetrahedralMesh* tetmesh ) const
{
	double R,r;
	
	R = tet_circ_r(tetmesh->m_storage_verts[Verts[0]].pos(), tetmesh->m_storage_verts[Verts[1]].pos(),
		tetmesh->m_storage_verts[Verts[2]].pos(), tetmesh->m_storage_verts[Verts[3]].pos());
	r = tet_in_r(tetmesh->m_storage_verts[Verts[0]].pos(), tetmesh->m_storage_verts[Verts[1]].pos(),
		tetmesh->m_storage_verts[Verts[2]].pos(), tetmesh->m_storage_verts[Verts[3]].pos());

	return 3*r/R;
	

}
double cwg::tet_circ_r( const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
											 const Cleaver::vec3& v2pos, const Cleaver::vec3& v3pos )
{
	Cleaver::vec3 v0 = v0pos;
	Cleaver::vec3 v1 = v1pos - v0pos;
	Cleaver::vec3 v2 = v2pos - v0pos;
	Cleaver::vec3 v3 = v3pos - v0pos;

	
	double	r = length((v1.dot(v1))*(v2.cross(v3)) + (v2.dot(v2))*(v3.cross(v1)) + (v3.dot(v3))*(v1.cross(v2)))
			/(2 * abs(dot(v1,cross(v2,v3))));
	
	

	return r;
}




double cwg::tet_in_r( const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
											 const Cleaver::vec3& v2pos, const Cleaver::vec3& v3pos )
{
	Cleaver::vec3 v0 = v0pos;
	Cleaver::vec3 v1 = v1pos - v0pos;
	Cleaver::vec3 v2 = v2pos - v0pos;
	Cleaver::vec3 v3 = v3pos - v0pos;

	double r = abs(dot(v1,cross(v2,v3))) / 
		(length(cross(v2,v3)) + length(cross(v3,v1)) + length(cross(v1,v2)) + 
		length(cross(v2,v3) + cross(v3,v1) + cross(v1,v2)) );
	return r;
}

double cwg::Tetrahedron::compute_delta_supervolume( const TetrahedralMesh* tetmesh ) const
{
	vector<Cleaver::vec3> verts_pos(4);
	for (int i=0; i<4; i++) 
		verts_pos[i] = tetmesh->m_storage_verts[ this->Verts[i] ].pos();
	return cwg::tet_delta_supervolume(verts_pos[0], verts_pos[1], verts_pos[2], verts_pos[3]);
}
double cwg::Tetrahedron::compute_delta_supervolume_sliding_result( const TetrahedralMesh* tetmesh ) const
{
	vector<Cleaver::vec3> verts_pos(4);
	for (int i=0; i<4; i++) 
		verts_pos[i] = tetmesh->ResultV[ this->Verts[i] ].pos();
	return cwg::tet_delta_supervolume(verts_pos[0], verts_pos[1], verts_pos[2], verts_pos[3]);
}
cwg::ordered_triple cwg::Tetrahedron::get_opposite_face( int cur_vid ) const
{
	int opfvs[3]; 
	int k = 0; 
	for (int i=0; i<4; i++)
	{
		if (this->Verts[i] == cur_vid)
			continue;

		opfvs[k++] = this->Verts[i];
	}
	return cwg::ordered_triple(opfvs[0], opfvs[1], opfvs[2]);
}

void cwg::shared_verts( Tetrahedron* ti, Tetrahedron* tj, std::vector<int>& titjverts )
{
	titjverts.clear();
	for (int i=0; i<4; i++)
	{
		for (int j=0; j<4; j++)
		{
			if (ti->Verts[i] == tj->Verts[j])
			{
				titjverts.push_back(ti->Verts[i]);
				break;
			}
		}
	}
}

std::ostream& cwg::operator<<( std::ostream& os, const cwg::Tetrahedron& obj )
{
	os << obj.ID << " ";
	for (int i=0; i<4; i++)
		os << obj.Verts[i] << " ";
	os << obj.m_label;
	return os;
}

std::istream& cwg::operator>>( std::istream& is, cwg::Tetrahedron& obj )
{
	is >> obj.ID >> obj.Verts[0] >> obj.Verts[1] >> obj.Verts[2] >> obj.Verts[3] >> obj.m_label;
	return is;
}

double cwg::tet_vol( const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
	const Cleaver::vec3& v2pos, const Cleaver::vec3& v3pos )
{
	return Cleaver::dot(v0pos - v3pos, cross(v1pos-v3pos, v2pos-v3pos)) / 6.0;
}

std::vector<double> cwg::tet_dihedral_angles( const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
	const Cleaver::vec3& v2pos, const Cleaver::vec3& v3pos )
{
	Cleaver::vec3 v01 = v1pos - v0pos;
	Cleaver::vec3 v02 = v2pos - v0pos;
	Cleaver::vec3 v03 = v3pos - v0pos;
	Cleaver::vec3 n012 = v01.cross(v02);
	Cleaver::vec3 v12 = v2pos - v1pos;
	Cleaver::vec3 v13 = v3pos - v1pos;

	vector<double> angles;
	const double precision_tmp = 1e-6;
	Cleaver::vec3 ns[4];
	if (v03.dot(n012) >= 0)
	{
		Cleaver::vec3 n021 = n012 * -1.0;
		Cleaver::vec3 n032 = v03.cross(v02);
		Cleaver::vec3 n013 = v01.cross(v03);
		Cleaver::vec3 n123 = v12.cross(v13);
		ns[0] = n021;
		ns[1] = n032;
		ns[2] = n013;
		ns[3] = n123;
	}
	else
	{
		Cleaver::vec3 n023 = v02.cross(v03);
		Cleaver::vec3 n031 = v03.cross(v01);
		Cleaver::vec3 n132 = v13.cross(v12);
		ns[0] = n012;
		ns[1] = n023;
		ns[2] = n031;
		ns[3] = n132;
	}

	for (int i=0; i<4; i++)
	{
		for (int j=i+1; j<4; j++)
		{
			double arccos_nij = ns[i].dot(ns[j]) / (Cleaver::length(ns[i])*Cleaver::length(ns[j]));
			if (abs(arccos_nij+1) < precision_tmp)
				angles.push_back(0.0);
			else if (abs(arccos_nij-1) < precision_tmp)
				angles.push_back(CV_PI);
			else
				angles.push_back(CV_PI - acos(arccos_nij));
		}
	}
	return angles;
}

double cwg::tet_delta_supervolume( const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
	const Cleaver::vec3& v2pos, const Cleaver::vec3& v3pos )
{
	hpan::vec3g<double> vec3g_v0pos(v0pos.x, v0pos.y, v0pos.z);
	hpan::vec3g<double> vec3g_v1pos(v1pos.x, v1pos.y, v1pos.z);
	hpan::vec3g<double> vec3g_v2pos(v2pos.x, v2pos.y, v2pos.z);
	hpan::vec3g<double> vec3g_v3pos(v3pos.x, v3pos.y, v3pos.z);

	double vol, svol0, svol1, svol2, svol3;
	svol0 = hpan::quadratic_integral_tet(vec3g_v0pos, vec3g_v0pos, vec3g_v1pos, vec3g_v2pos, vec3g_v3pos, vol);
	svol1 = hpan::quadratic_integral_tet(vec3g_v1pos, vec3g_v0pos, vec3g_v1pos, vec3g_v2pos, vec3g_v3pos, vol);
	svol2 = hpan::quadratic_integral_tet(vec3g_v2pos, vec3g_v0pos, vec3g_v1pos, vec3g_v2pos, vec3g_v3pos, vol);
	svol3 = hpan::quadratic_integral_tet(vec3g_v3pos, vec3g_v0pos, vec3g_v1pos, vec3g_v2pos, vec3g_v3pos, vol);
	return abs( svol0+svol1+svol2+svol3 ) / 6.0;
}

