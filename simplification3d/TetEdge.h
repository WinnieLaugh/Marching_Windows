#ifndef _TET_EDGE_H_
#define _TET_EDGE_H_

#include "TetVertex.h"
#include <fstream>
#include <set>
namespace cwg
{
	class TetEdgeAroundInfo
	{
	public:
		inline void						clear() { edge_tets.clear(); v0_private_tets.clear(); v1_private_tets.clear();
												v0_private_overlap_tets.clear(); v1_private_overlap_tets.clear(); 
												v0_private_nonoverlap_tets.clear(); v1_private_nonoverlap_tets.clear(); 
												affected_verts.clear(); };

		std::vector<int>				edge_tets;
		std::vector<int>				v0_private_tets;
		std::vector<int>				v1_private_tets;
		std::vector<int>				v0_private_overlap_tets;
		std::vector<int>				v0_private_nonoverlap_tets;
		std::vector<int>				v1_private_overlap_tets;
		std::vector<int>				v1_private_nonoverlap_tets;

		std::set<int>					affected_verts;
	};

	class TetEdge
	{
	public:
		TetEdge() : ID(ms_nextid++), onBoundary(false), able_to_collapse(false), Cost(-1.) {}
		TetEdge(int vi, int vj) : ID(ms_nextid++), onBoundary(false), able_to_collapse(false) { Verts[0] = vi; Verts[1] = vj; }
		~TetEdge() {}

		int								find_vert(int vid) const;
		void							ComputeContraction(TetrahedralMesh* tetmesh);
		void							ComputeCost_Q_s(TetrahedralMesh* tetmesh);
		void							ComputeCost_Q_d(TetrahedralMesh* tetmesh);
		//===================================
		void							interface_weight(TetrahedralMesh* tetmesh, std::set<int>& interV);
		
		//=====================================

		inline void						set_verts(int vi, int vj) { Verts[0] = vi; Verts[1] = vj; }

		inline int						id() const { return ID; }
		inline double					cost() const { return Cost; }

		int								ID;
		int								Verts[2];
		Cleaver::vec3					ContractPos;
		double							Cost;
		bool							onBoundary;
		bool							able_to_collapse;
		bool							debug_bool;

		// for surface improvement
		//std::vector<int>				tri; // the triangles shared this edge
		int								tri[2];

		void							compute_around_info(const TetrahedralMesh* tetmesh);
		void							debug_obj(const TetrahedralMesh* tetmesh, std::string filename, int tet_index);
		void							debug_edge(const TetrahedralMesh* tetmesh, std::string filename, int vert_0, int vert_1);

		bool							validate_folderover(TetrahedralMesh* tetmesh);
		/////////////////////////////////// ********************** DEBUG *********************** ///////////////////
		bool							check_labels_debug(TetrahedralMesh* tetmesh, int label_id);
		bool							check_tet_debug(TetrahedralMesh* tetmesh, int tet_id);
		void							set_debug_bool(bool to_debug_bool);
		/////////////////////////////////// ********************** DEBUG *********************** ///////////////////

		inline const TetEdgeAroundInfo&	get_around_info() const { return ms_edge_around_info; }
		inline static void				reset_nextid() { ms_nextid = 0; }

	private:
		void							qem_pos(const cv::Matx33d& Q3x3, const cv::Vec3d& q3x1, TetVertex* v0, TetVertex* v1);
		int								get_order( const TetrahedralMesh* tetmesh );
		int								get_edge_order( const TetrahedralMesh* tetmesh, int v0, int v1);
		bool							check_dual_edges( const TetrahedralMesh* tetmesh);
		bool							check_dual_faces( const TetrahedralMesh* tetmesh, int v0, int v1);
		bool							check_volume_dual_faces( const TetrahedralMesh* tetmesh, int v0, int v1);
		inline void						qem_pos(double diagval, const Cleaver::vec3& q3x1) { ContractPos = q3x1 * (-1.0) / diagval; }
		bool							check_edge_pos(const TetrahedralMesh* tetmesh, double x_0, double y_0, double z_0, double x_1, double y_1, double z_1);


		static int						ms_nextid;

		static TetEdgeAroundInfo		ms_edge_around_info;
	};

}
#endif