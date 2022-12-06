#ifndef _TETRAHEDRON_H_
#define _TETRAHEDRON_H_

#ifdef SIMPLIFICATION3D_EXPORTS
#define SIMPLIFICATION3D_API __declspec(dllexport)
#else
#define SIMPLIFICATION3D_API __declspec(dllimport)
#endif

#include "TetVertex.h"
#include "TetEdge.h"
#include "TetTriangle.h"
#include <Data.h>

namespace cwg
{
	class TetrahedralMesh;

	class SIMPLIFICATION3D_API Tetrahedron : public Data
	{
	public:
		Tetrahedron() : ID(ms_nextid++){Data();}
		Tetrahedron(int id, int vids[4], int lb) : ID(id), m_label(lb) {Data(); for (int i=0; i<4; i++) Verts[i] = vids[i]; }

		~Tetrahedron() {}

		int								ID;
		bool							Verts_last_slide[4];	//if the v is in last slide moving part
		double							variance(const TetrahedralMesh* tetmesh, double avg) const;
		double							avg_edge_length_sliding_result(const TetrahedralMesh* tetmesh) const;

		void							correlate_verts(TetrahedralMesh* tetmesh);
		void							correlate_verts_sliding_result(TetrahedralMesh* tetmesh);
		int								find_vert(int vid) const;
		bool							inside_tet_check(float x, float y, float z, std::vector<TetVertex>	*vertex_vector);

		ordered_triple					get_opposite_face(int cur_vid) const;

		double							compute_volume(const TetrahedralMesh* tetmesh) const;
		double compute_volume_sliding_result(const TetrahedralMesh* tetmesh) const;
		std::vector<double>				compute_dihedral_angle(const TetrahedralMesh* tetmesh) const;
		double							minimal_dihedral_angle(const TetrahedralMesh* tetmesh) const;
		double							compute_delta_supervolume(const TetrahedralMesh* tetmesh) const;
		double							compute_delta_supervolume_sliding_result(const TetrahedralMesh* tetmesh) const;
		void							compute_sphere_ratio(TetrahedralMesh* tetmesh) const;
		double							compute_tet_ratio(TetrahedralMesh* tetmesh) const;
		void							compute_tri_ratio(TetrahedralMesh* tetmesh, std::vector<double>& surf_ratio) const;
		void							compute_tri_ratio_cgal(TetrahedralMesh* tetmesh, std::vector<double>& surf_ratio) const;

		inline void						set_label(int lb) { m_label = lb; }
		inline int						get_label() const { return m_label; }

		inline void						set_ratio(double r) { r_ratio = r; }
		inline double					get_ratio() const { return r_ratio; } 

		void							read(std::ifstream& ifile);
		void							write(std::ofstream& ofile) const;
		inline static void				reset_nextid() { ms_nextid = 0; }

		friend std::ostream&			operator<<(std::ostream& os, const Tetrahedron& obj);
		friend std::istream&			operator>>(std::istream& is, Tetrahedron& obj);

		inline void						set_min_angle(float min_angle) { min_dihedral_angle = min_angle; }
		inline float					get_min_angle() const { return min_dihedral_angle; } 
		void							set_Verts(int verts[4]){ Verts[0] = verts[0]; Verts[1] = verts[1]; Verts[2] = verts[2]; Verts[3] = verts[3];}
		void							invalidate(TetrahedralMesh* tetmesh);
	private:
		int								m_label;
		float							min_dihedral_angle;
		double							r_ratio;
		
	private:
		static int						ms_nextid;
	};

	double								tet_vol(const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
											const Cleaver::vec3& v2pos, const Cleaver::vec3& v3pos);
	std::vector<double>					tet_dihedral_angles(const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
											const Cleaver::vec3& v2pos, const Cleaver::vec3& v3pos);
	double								tet_delta_supervolume(const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
											const Cleaver::vec3& v2pos, const Cleaver::vec3& v3pos);

	double								tet_circ_r(const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
											const Cleaver::vec3& v2pos, const Cleaver::vec3& v3pos);	//radiu of circumsphere
	double								tet_in_r(const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
											const Cleaver::vec3& v2pos, const Cleaver::vec3& v3pos);	//radiu of insphere

	Cleaver::vec3&						tet_circ_c(const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
											const Cleaver::vec3& v2pos, const Cleaver::vec3& v3pos);


	double								tri_circ_r(const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
											const Cleaver::vec3& v2pos);	//radiu of circumcircle
	double								tri_in_r(const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
											const Cleaver::vec3& v2pos);	//radiu of incircle
	double								tri_area(const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
											const Cleaver::vec3& v2pos);	//area of triangle

	Cleaver::vec3&						tri_circ_c(const Cleaver::vec3& v0pos, const Cleaver::vec3& v1pos, 
											const Cleaver::vec3& v2pos);	//center of circumcircle

	void								shared_verts(Tetrahedron* ti, Tetrahedron* tj, std::vector<int>& titjverts);
	std::ostream&						operator<<(std::ostream& os, const Tetrahedron& obj);
	std::istream&						operator>>(std::istream& is, Tetrahedron& obj);
}

#endif